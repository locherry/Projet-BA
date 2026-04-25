from __future__ import annotations

from typing import Dict

from .Materials import Concrete, Steel
from .TSection import TSection
from .ReinforcementLayout import ReinforcementLayout
import numpy as np
from scipy.optimize import brentq


class FlexuralDesign:
    """
    EC2 flexural design and SLS checks for a T-section beam.

    All calculations use a single consistent effective depth:
      • At construction d is a first estimate (e.g. 0.9·h).
      • Once a ReinforcementLayout is attached, call update_d() to replace d
        with d_reel so that ULS sizing, SLS stresses, and crack control all
        operate on the same lever arm.
    """

    def __init__(
        self,
        section: TSection,
        concrete: Concrete,
        steel: Steel,
        d: float,
    ) -> None:
        self.section = section
        self.concrete = concrete
        self.steel = steel
        self.d = d                          # updated via update_d() after layout is set
        self.layout: ReinforcementLayout | None = None

    # ---------------------------------------------------------------------- #
    #  Layout attachment                                                        #
    # ---------------------------------------------------------------------- #

    def set_layout(self, layout: ReinforcementLayout) -> None:
        """Attach a ReinforcementLayout and immediately update d to d_reel."""
        self.layout = layout
        self.update_d()

    def update_d(self) -> None:
        """Replace self.d with the area-weighted d_reel from the current layout."""
        if self.layout is None:
            raise ValueError("No layout attached — call set_layout() first.")
        self.d = self.layout.compute_d_reel()

    # ---------------------------------------------------------------------- #
    #  ULS — neutral axis and reinforcement area                               #
    # ---------------------------------------------------------------------- #

    def neutral_axis(self, M_Ed: float) -> Dict[str, float]:
        """
        ULS neutral axis depth using the EC2 rectangular stress block.
        Currently implemented for neutral axis in the flange only.
        Raises NotImplementedError with a clear message when the NA falls in the web.
        """
        b_eff = self.section.b_eff
        h_f = self.section.h_f
        f_cd = self.concrete.f_cd
        d = self.d

        # Moment capacity with full flange in compression
        M_flange = b_eff * h_f * f_cd * (d - h_f / 2)

        if M_Ed > M_flange:
            raise NotImplementedError(
                f"M_Ed ({M_Ed/1e6:.1f} MN·m) exceeds flange capacity "
                f"({M_flange/1e6:.1f} MN·m) — neutral axis in web: "
                "implement composite T-section ULS before proceeding."
            )

        mu = M_Ed / (b_eff * d ** 2 * f_cd)
        alpha = (1 - np.sqrt(1 - 2 * mu)) / 0.8
        x_na = alpha * d
        return {"x_na": x_na, "alpha": alpha, "mu": mu}

    def reinforcement_area(self, M_Ed: float) -> float:
        """
        Required tensile reinforcement area for ULS moment M_Ed.
        Uses self.d (updated to d_reel once layout is attached).
        """
        res = self.neutral_axis(M_Ed)
        alpha = res["alpha"]

        epsilon_c = 3.5e-3
        epsilon_s = epsilon_c * (1 - alpha) / alpha
        sigma_s = self.steel.stress(epsilon_s)

        F_c = 0.8 * alpha * self.section.b_eff * self.d * self.concrete.f_cd
        A_s = F_c / sigma_s
        return A_s

    # ---------------------------------------------------------------------- #
    #  Creep / modular ratio helpers                                           #
    # ---------------------------------------------------------------------- #

    def _creep_coefficient(self, t0: float = 28.0) -> float:
        """
        φ(∞, t0) from EC2 Annex B.
        Uses RH = 50 %, h0 computed from b_w and h_tot.
        t0: age at loading in days.
        """
        f_cm_MPa = self.concrete.f_cm * 1e-6        # Pa → MPa

        a1 = (35 / f_cm_MPa) ** 0.7
        a2 = (35 / f_cm_MPa) ** 0.2

        # Notional size h0 = 2·Ac / u  (mm) — using web rectangle as reference
        h0 = (
            (self.section.b_w * self.section.h_tot)
            / (self.section.b_w + self.section.h_tot)
            * 1e3          # → mm
        )

        RH = 0.50          # 50 % relative humidity
        phi_RH = (1 + (1 - RH) / (0.1 * h0 ** (1 / 3)) * a1) * a2
        beta_fcm = 16.8 / np.sqrt(f_cm_MPa)
        beta_t0 = 1.0 / (0.1 + t0 ** 0.2)
        
        return phi_RH * beta_fcm * beta_t0

    def _phi_effective(self, phi_inf_t0: float, M_Eqp: float, M_Ed: float) -> float:
        """Effective creep coefficient for SLS: φ_ef = φ(∞,t0)·M_Eqp/M_Ed."""
        if M_Ed == 0:
            raise ValueError("M_Ed must be non-zero to compute φ_ef")
        return phi_inf_t0 * (M_Eqp / M_Ed)

    def _cracked_na_sls(self, A_s: float, n: float) -> float:
        """
        Cracked (homogeneous) neutral axis depth x for T-section under SLS.

        Equilibrium (first moment of area about NA = 0):
          b_w·x²/2 + (b_eff − b_w)·h_f·(x − h_f/2) − n·A_s·(d − x) = 0

        Solved analytically as a quadratic:
          a = b_w/2
          b = (b_eff − b_w)·h_f + n·A_s
          c = −(b_eff − b_w)·h_f·(h_f/2) − n·A_s·d
        """
        b_w = self.section.b_w
        b_eff = self.section.b_eff
        h_f = self.section.h_f
        d = self.d

        a = b_w / 2
        b = (b_eff - b_w) * h_f + n * A_s
        c = -(b_eff - b_w) * h_f * (h_f / 2) - n * A_s * d

        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            raise ValueError("No real solution for cracked neutral axis.")
        x = (-b + np.sqrt(discriminant)) / (2 * a)
        return x

    def _cracked_inertia(self, x: float, A_s: float, n: float) -> float:
        """
        Cracked second moment of area I_hr about the cracked NA.

        • x ≤ h_f  → NA in flange: rectangular section b_eff
        • x  > h_f → NA in web:    composite T geometry
        """
        b_w = self.section.b_w
        b_eff = self.section.b_eff
        h_f = self.section.h_f
        d = self.d

        if x <= h_f:
            I_hr = b_eff * x ** 3 / 3 + n * A_s * (d - x) ** 2
        else:
            I_hr = (
                b_w * x ** 3 / 3
                + (b_eff - b_w) * h_f ** 3 / 12
                + (b_eff - b_w) * h_f * (x - h_f / 2) ** 2
                + n * A_s * (d - x) ** 2
            )
        return I_hr

    # ---------------------------------------------------------------------- #
    #  SLS — stress limitation  (EC2 §7.2)                                    #
    # ---------------------------------------------------------------------- #

    def check_stress_limitation(
        self,
        M_ELS_car: float,
        M_ELS_eqp: float,
        A_s: float,
        t0: float = 28,
    ) -> dict:
        """
        SLS stress limitation check per EC2 §7.2.

        Limits
        ------
        Concrete compression (quasi-permanent):  σ_c ≤ 0.45·f_ck  (§7.2(3))
        Steel tension        (characteristic):   σ_s ≤ 0.80·f_yk  (§7.2(5))

        Parameters
        ----------
        M_ELS_car  : ULS design moment [N·m]  — used to derive φ_ef
        M_ELS_eqp : quasi-permanent SLS moment [N·m] — used for stress calculation
        A_s   : reinforcement area [m²]
        t0    : age at loading [days]
        """
        if self.layout is None:
            raise ValueError("Layout must be set via set_layout() before SLS checks.")

        # ---- creep and effective modular ratio ----
        phi_inf_t0 = self._creep_coefficient(t0)
        phi_ef = self._phi_effective(phi_inf_t0, M_ELS_eqp, M_ELS_car)
        E_c_eff = self.concrete.E_cm / (1.0 + phi_ef)
        n = self.steel.Es / E_c_eff          # αe with creep

        # ---- cracked NA and inertia (single shared method) ----
        x = self._cracked_na_sls(A_s, n)
        I_hr = self._cracked_inertia(x, A_s, n)

        # ---- stresses under quasi-permanent moment ----
        sigma_c = M_ELS_eqp * x / I_hr                    # concrete top fibre
        sigma_s = n * M_ELS_eqp * (self.d - x) / I_hr    # steel

        # ---- EC2 limits ----
        limit_compression = 0.45 * self.concrete.f_ck   # §7.2(3) quasi-permanent
        limit_steel = 0.80 * self.steel.f_yk             # §7.2(5) characteristic

        return {
            "phi_inf_t0": phi_inf_t0,
            "phi_ef": phi_ef,
            "E_cm": self.concrete.E_cm,
            "E_c_eff": E_c_eff,
            "n": n,
            "x_SLS": x,
            "I_hr": I_hr,
            "sigma_c": sigma_c,
            "sigma_s": sigma_s,
            "limit_compression": limit_compression,
            "limit_steel": limit_steel,
            "ok_compression": sigma_c <= limit_compression,
            "ok_steel": sigma_s <= limit_steel,
        }

    # ---------------------------------------------------------------------- #
    #  SLS — crack control  (EC2 §7.3)                                        #
    # ---------------------------------------------------------------------- #

    def crack_control(
        self,
        M_ELS_car: float,
        M_ELS_eqp: float,
        A_s: float,
        w_max: float = 0.4e-3,
    ) -> dict:
        """
        Crack control check per EC2 §7.3.

        Parameters
        ----------
        M_ELS_car  : ULS design moment [N·m]
        M_ELS_eqp : quasi-permanent SLS moment [N·m]
        A_s   : provided reinforcement area [m²]
        w_max : maximum crack width [m] from EC2 Table 7.1N for the exposure class.
                Default 0.4 mm (XC1, quasi-permanent).  Pass 0.3e-3 or 0.2e-3
                for more aggressive exposures.
        """
        if self.layout is None:
            raise ValueError("Layout must be set via set_layout() before SLS checks.")

        b_w = self.section.b_w
        b_eff = self.section.b_eff
        h_f = self.section.h_f
        h = self.section.h_tot
        d = self.d
        n = self.steel.Es / self.concrete.E_cm      # αe without creep for crack control
        f_ct_eff = self.concrete.f_ctm
        k_t = 0.4                                    # long-term loading (EC2 §7.3.4)

        # ------------------------------------------------------------------ #
        #  1. Non-linear stress distribution factor k  (EC2 §7.3.2)          #
        # ------------------------------------------------------------------ #
        def k_factor(x_na: float) -> float:
            """
            k accounts for non-uniform self-equilibrating stresses.
            h or b_eff is the relevant dimension depending on which zone is in tension.
            """
            dim = b_eff if x_na < h_f else h
            if dim < 0.300:
                return 1.0
            elif dim >= 0.800:
                return 0.65
            else:
                return 1.21 - 0.70 * dim

        # ------------------------------------------------------------------ #
        #  2. Minimum reinforcement  (EC2 §7.3.2)                             #
        # ------------------------------------------------------------------ #
        k_c = 0.4                  # pure bending
        x_na_uls = self.neutral_axis(M_ELS_car)["x_na"]
        A_ct = b_w * (h - x_na_uls)          # tensile concrete area before cracking
        As_min = k_c * k_factor(x_na_uls) * f_ct_eff * A_ct / self.steel.f_yk

        # ------------------------------------------------------------------ #
        #  3. Cracked NA under SLS (αe = E_s/E_cm, no creep for crack)       #
        # ------------------------------------------------------------------ #
        x = self._cracked_na_sls(A_s, n)
        I_hr = self._cracked_inertia(x, A_s, n)

        # ------------------------------------------------------------------ #
        #  4. Steel stress under quasi-permanent combination                  #
        # ------------------------------------------------------------------ #
        sigma_s = n * M_ELS_eqp * (d - x) / I_hr

        # ------------------------------------------------------------------ #
        #  5. Effective concrete area and reinforcement ratio (EC2 §7.3.4)   #
        #     h_c,ef = min{ 2.5(h−d) ; (h−x)/3 ; h/2 }                      #
        # ------------------------------------------------------------------ #
        h_c_ef = min(2.5 * (h - d), (h - x) / 3, h / 2)
        A_c_eff = b_w * h_c_ef
        rho_p_eff = A_s / A_c_eff

        # ------------------------------------------------------------------ #
        #  6. Mean strain difference ε_sm − ε_cm  (EC2 §7.3.4(2))           #
        #     ≥ 0.6·σ_s/E_s                                                  #
        # ------------------------------------------------------------------ #
        eps_diff = (
            sigma_s - k_t * (f_ct_eff / rho_p_eff) * (1.0 + n * rho_p_eff)
        ) / self.steel.Es
        eps_diff = max(eps_diff, 0.6 * sigma_s / self.steel.Es)

        # ------------------------------------------------------------------ #
        #  7. Maximum crack spacing s_r,max  (EC2 §7.3.4(3))                 #
        #     k1=0.8 (HA bars), k2=0.5 (bending), k4=0.425                  #
        # ------------------------------------------------------------------ #
        k1, k2, k4 = 0.8, 0.5, 0.425
        c_nom_mm = self.layout.c_nom * 1e3
        k3 = 3.4 if c_nom_mm <= 25 else 3.4 * (25 / c_nom_mm) ** (2 / 3)

        phi_eq = self.layout.phi_eq()
        s_r_max = k3 * self.layout.c_nom + k1 * k2 * k4 * phi_eq / rho_p_eff

        # ------------------------------------------------------------------ #
        #  8. Crack width  w_k = s_r,max · (ε_sm − ε_cm)  (EC2 §7.3.4(1))  #
        # ------------------------------------------------------------------ #
        w_k = s_r_max * eps_diff

        return {
            # Min reinforcement
            "As_min": As_min,
            "ok_As_min": A_s >= As_min,
            # SLS section state
            "x_SLS": x,
            "I_hr": I_hr,
            "sigma_s": sigma_s,
            # Effective area
            "h_c_ef": h_c_ef,
            "A_c_eff": A_c_eff,
            "rho_p_eff": rho_p_eff,
            # Strain difference
            "eps_diff": eps_diff,
            # Crack spacing and width
            "phi_eq": phi_eq,
            "s_r_max": s_r_max,
            "w_k": w_k,
            "w_max": w_max,
            "ok_wk": w_k <= w_max,
        }