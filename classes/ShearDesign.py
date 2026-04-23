from typing import Tuple, Optional, cast
from .TSection import TSection
import numpy as np


class ShearDesign:
    """
    Dimensionne les armatures d'effort tranchant (cadres/épingles)
    selon l'analogie du treillis de Ritter-Mörsch (EC2 §6.2).

    Hypothèses
    ----------
    * Cadres verticaux (α = 90°)
    * Angle de bielle θ tel que cot(θ) ∈ [1 ; 2.5]  (EC2 §6.2.3(2))
    * Bras de levier interne z = 0.9 d
    """

    # EC2 §6.2.3(2) : 1 ≤ cot θ ≤ 2.5
    COT_THETA_MIN = 1.0
    COT_THETA_MAX = 2.5

    def __init__(
        self,
        section: TSection,
        v_ed: np.ndarray,
        x: np.ndarray,
        d: float,
        f_ywd: float,
        f_cd: float,
        phi_t: float = 8e-3,
    ):
        """
        Parameters
        ----------
        v_ed  : efforts tranchants ELU V_Ed(x) [N]
        x     : abscisses [m]  (même longueur que v_ed)
        d     : hauteur utile [m]
        f_ywd : résistance de calcul acier transversal [Pa]
        f_cd  : résistance de calcul béton [Pa]
        phi_t : diamètre des cadres [m]
        """
        self.bw = section.b_w
        self.section = section
        self.v_ed  = np.asarray(v_ed,  dtype=float)
        self.x     = np.asarray(x,     dtype=float)
        self.d     = float(d)
        self.f_ywd = float(f_ywd)
        self.f_cd  = float(f_cd)
        self.phi_t = float(phi_t)

        self.area_t_per_leg = np.pi * (phi_t / 2.0) ** 2  # section d'un brin [m²]
        self.z = 0.9 * d                                   # bras de levier interne [m]

    # ------------------------------------------------------------------ #
    #  Helpers
    # ------------------------------------------------------------------ #

    @property
    def bw(self) -> Optional[float]:
        """Largeur de l'âme — doit être définie par sous-classe ou setattr."""
        return getattr(self, "_bw", None)

    @bw.setter
    def bw(self, value: float):
        self._bw = float(value)

    def _check_bw(self):
        if self.bw is None:
            raise AttributeError(
                "bw (largeur d'âme) doit être renseigné avant d'appeler "
                "cette méthode : shear.bw = section.b_w"
            )

    # ------------------------------------------------------------------ #
    #  EC2 §6.2.2 — Résistance sans armatures transversales V_Rd,c
    # ------------------------------------------------------------------ #

    def v_rd_c(self, A_sl: float, N_Ed: float = 0.0) -> np.ndarray:
        """
        Résistance au cisaillement sans armatures transversales V_Rd,c [N].

        EC2 eq. (6.2a) / (6.2b) — valeur scalaire ou tableau selon x.

        Parameters
        ----------
        A_sl : section des armatures longitudinales tendues [m²]
        N_Ed : effort normal de calcul [N]  (positif = compression)
        """
        self._check_bw()
        bw = cast(float, self.bw)
        d  = self.d

        # EC2 eq. (6.2) — facteur taille
        k = np.minimum(1.0 + np.sqrt(200e-3 / d), 2.0)

        # Taux de ferraillage longitudinal (borné à 0.02)
        rho_l = np.minimum(A_sl / (bw * d), 0.02)

        # f_ck en MPa pour la formule EC2
        f_ck_MPa = self.f_cd * 1.5 / 1e6          # on remonte à f_ck [MPa]

        # Coefficient C_Rd,c
        C_Rd_c = 0.18 / 1.5                        # = 0.12

        # v_min  EC2 eq. (6.3N)
        v_min = 0.035 * k ** 1.5 * np.sqrt(f_ck_MPa)   # [MPa^0.5 * dim'less]

        # Contrainte normale [MPa]
        
        # FIX: use gross area A_c of the T-section, not a derived area
        A_c = self.section.b_eff * self.section.h_f + self.section.b_w * self.section.h_w
        sigma_cp = min(N_Ed / A_c, 0.2 * f_ck_MPa * 1e6) / 1e6  # [MPa]

        # sigma_cp = min(N_Ed / (bw * (d / 0.9)), 0.2 * f_ck_MPa * 1e6) / 1e6

        # V_Rd,c eq. (6.2a)
        v_rdc_a = (C_Rd_c * k * (100 * rho_l * f_ck_MPa) ** (1 / 3)
                   + 0.15 * sigma_cp) * bw * d * 1e6   # [N]

        # V_Rd,c eq. (6.2b) — valeur minimale
        v_rdc_b = (v_min + 0.15 * sigma_cp) * bw * d * 1e6

        return np.maximum(v_rdc_a, v_rdc_b)

    # ------------------------------------------------------------------ #
    #  EC2 §6.2.3 — Angle de bielle optimal θ
    # ------------------------------------------------------------------ #

    def cot_theta_optimal(self) -> float:
        """
        Angle de bielle θ optimal selon EC2 §6.2.3(2).

        On choisit cot θ = 2.5 (θ ≈ 21.8°) pour minimiser la quantité
        d'armatures transversales, sous réserve de la résistance de la bielle.
        La vérification de V_Rd,max est faite dans `v_rd_max`.
        """
        return self.COT_THETA_MAX   # valeur la plus économique

    def cot_theta_required(self) -> float:
        """
        Valeur minimale de cot θ imposée par V_Ed max (EC2 eq. 6.9).
        Si V_Ed > V_Rd,max(cot_θ=1), le dimensionnement est impossible.
        """
        self._check_bw()
        V_Ed_max = np.max(np.abs(self.v_ed))
        bw = cast(float, self.bw)
        z  = self.z

        # ν = 0.6*(1 - f_ck/250)  avec f_ck en MPa
        f_ck_MPa = self.f_cd * 1.5 / 1e6
        nu = 0.6 * (1.0 - f_ck_MPa / 250.0)

        # V_Rd,max(cot θ) = ν·f_cd·bw·z / (cot θ + tan θ)
        # The denominator (cot θ + 1/cot θ) has a minimum at cot θ = 1 (θ = 45°),
        # so V_Rd,max is highest at cot θ = 1 and decreases toward cot θ = 2.5.
        # Scanning from MIN to MAX finds the smallest cot θ that still satisfies
        # V_Rd,max ≥ V_Ed,max — i.e. the least-inclined strut that is still safe.
        cot_vals = np.linspace(self.COT_THETA_MIN, self.COT_THETA_MAX, 2000)
        tan_vals = 1.0 / cot_vals
        v_rdmax_vals = nu * self.f_cd * bw * z / (cot_vals + tan_vals)

        # Indice où V_Rd,max ≥ V_Ed_max
        ok = np.where(v_rdmax_vals >= V_Ed_max)[0]
        if len(ok) == 0:
            raise ValueError(
                f"V_Ed,max = {V_Ed_max/1e3:.1f} kN dépasse V_Rd,max même "
                f"avec cot θ = {self.COT_THETA_MIN:.1f}. "
                "Augmenter bw ou f_ck."
            )
        # cot θ minimal qui satisfait la condition
        return float(cot_vals[ok[0]])

    # ------------------------------------------------------------------ #
    #  EC2 §6.2.3 eq. (6.9) — Résistance maximale des bielles V_Rd,max
    # ------------------------------------------------------------------ #

    def v_rd_max(self, cot_theta: float | None = None) -> np.ndarray:
        """
        Résistance maximale de la bielle comprimée V_Rd,max [N].

        EC2 eq. (6.9) :
            V_Rd,max = α_cw · bw · z · ν₁ · f_cd / (cot θ + tan θ)

        Parameters
        ----------
        cot_theta : si None, utilise cot_theta_optimal()
        """
        self._check_bw()
        if cot_theta is None:
            cot_theta = self.cot_theta_optimal()

        f_ck_MPa = self.f_cd * 1.5 / 1e6
        nu1  = 0.6 * (1.0 - f_ck_MPa / 250.0)   # EC2 §6.2.3(3)
        alpha_cw = 1.0                             # pas de précontrainte

        bw = cast(float, self.bw)
        v_rdmax = (alpha_cw * bw * self.z * nu1 * self.f_cd
                   / (cot_theta + 1.0 / cot_theta))

        # Scalar → broadcast to x shape
        return np.full_like(self.x, v_rdmax)

    # ------------------------------------------------------------------ #
    #  EC2 §6.2.3 eq. (6.8) — Densité d'armatures transversales Asw / st
    # ------------------------------------------------------------------ #

    def asw_per_spacing(self, cot_theta: float | None = None) -> np.ndarray:
        """
        Densité d'armatures transversales requise A_sw/s_t [m²/m] le long de x.

        EC2 eq. (6.8) :
            A_sw / s_t = V_Ed / (z · f_ywd · cot θ)

        Parameters
        ----------
        cot_theta : si None, utilise cot_theta_optimal()
        """
        if cot_theta is None:
            cot_theta = self.cot_theta_optimal()

        return np.abs(self.v_ed) / (self.z * self.f_ywd * cot_theta)

    # ------------------------------------------------------------------ #
    #  Calcul de l'espacement des cadres
    # ------------------------------------------------------------------ #

    def stirrup_spacing(
        self,
        n_legs: int = 2,
        cot_theta: float | None = None,
    ) -> np.ndarray:
        """
        Espacement des cadres s_t [m] pour n_legs brins par cadre.

        s_t = (n_legs · A_t,leg) / (A_sw / s_t)

        Parameters
        ----------
        n_legs    : nombre de brins par cadre (2 pour un cadre simple)
        cot_theta : angle de bielle (None → optimal)
        """
        asw_s = self.asw_per_spacing(cot_theta)
        # Eviter division par zéro aux zones de V_Ed ≈ 0
        with np.errstate(divide="ignore", invalid="ignore"):
            spacing = np.where(
                asw_s > 0,
                n_legs * self.area_t_per_leg / asw_s,
                np.inf,
            )
        return spacing

    # ------------------------------------------------------------------ #
    #  EC2 §9.2.2 — Armature minimale transversale
    # ------------------------------------------------------------------ #

    def asw_min_per_spacing(self) -> float:
        """
        Densité minimale d'armatures transversales A_sw,min / s_t [m²/m].

        EC2 §9.2.2(5) eq. (9.4) :
            ρ_w,min = 0.08 · √f_ck / f_yk
            A_sw,min / s_t = ρ_w,min · bw · sin α   (α=90° → sin α=1)
        """
        self._check_bw()
        f_ck_MPa = self.f_cd * 1.5 / 1e6
        # f_ywd is assumed to be f_yk / γ_s (= f_yk / 1.15); we recover f_yk here.
        # If f_ywd was passed as the characteristic value, this is wrong — assert guards it.
        assert self.f_ywd < 600e6, (
            "f_ywd appears to be a characteristic value, not a design value "
            "(expected f_yk/1.15). Pass f_yd = f_yk/1.15."
        )
        f_yk_MPa = self.f_ywd * 1.15 / 1e6
        rho_w_min = 0.08 * np.sqrt(f_ck_MPa) / f_yk_MPa
        bw = cast(float, self.bw)
        return rho_w_min * bw

    # ------------------------------------------------------------------ #
    #  EC2 §9.2.2 — Espacement maximal des cadres
    # ------------------------------------------------------------------ #

    def max_stirrup_spacing(self, cot_theta: float | None = None) -> float:
        """
        Espacement maximal des cadres s_t,max [m].

        EC2 §9.2.2(6) :
            s_t,max = 0.75 · d · (1 + cot α)   avec α=90° → cot α=0
        """
        if cot_theta is None:
            cot_theta = self.cot_theta_optimal()
        return 0.75 * self.d   # cadres verticaux

    # ------------------------------------------------------------------ #
    #  Normalisation de l'espacement à une série constructive
    # ------------------------------------------------------------------ #

    @staticmethod
    def round_down_spacing(s: float, series_mm: list[int] | None = None) -> float:
        """
        Arrondit s à la valeur inférieure la plus proche d'une série
        constructive standard [mm].  Retourne la valeur en [m].
        """
        if series_mm is None:
            series_mm = [50, 75, 100, 125, 150, 175, 200, 225, 250, 300]
        s_mm = s * 1e3
        valid = [v for v in series_mm if v <= s_mm]
        if not valid:
            return series_mm[0] * 1e-3
        return max(valid) * 1e-3

    # ------------------------------------------------------------------ #
    #  Résumé du dimensionnement
    # ------------------------------------------------------------------ #

    def design_summary(
        self,
        n_legs: int = 2,
        cot_theta: float | None = None,
    ) -> dict:
        """
        Retourne un dictionnaire complet du dimensionnement en cisaillement.

        Returns
        -------
        dict avec les clés :
            cot_theta, theta_deg,
            v_ed_max, v_rd_max,
            asw_s_max  [m²/m],
            s_required_min [m],
            s_max [m],
            s_adopted [m],
            asw_s_provided [m²/m],
            ok_strut, ok_spacing
        """
        if cot_theta is None:
            cot_theta = self.cot_theta_optimal()

        theta_deg = np.degrees(np.arctan(1.0 / cot_theta))
        v_ed_max  = float(np.max(np.abs(self.v_ed)))
        vrdmax    = float(self.v_rd_max(cot_theta)[0])    # scalaire

        asw_s     = self.asw_per_spacing(cot_theta)
        asw_s_max = float(np.max(asw_s))

        s_req     = float(np.min(self.stirrup_spacing(n_legs, cot_theta)
                                 [np.isfinite(self.stirrup_spacing(n_legs, cot_theta))]))
        s_max     = self.max_stirrup_spacing(cot_theta)

        # Espacement adopté = min(s_req arrondi, s_max)
        s_adopted = min(self.round_down_spacing(s_req), s_max)

        asw_s_provided = n_legs * self.area_t_per_leg / s_adopted

        return dict(
            cot_theta       = cot_theta,
            theta_deg       = theta_deg,
            v_ed_max        = v_ed_max,
            v_rd_max        = vrdmax,
            asw_s_max       = asw_s_max,
            s_required_min  = s_req,
            s_max           = s_max,
            s_adopted       = s_adopted,
            asw_s_provided  = asw_s_provided,
            ok_strut        = v_ed_max <= vrdmax,
            ok_spacing      = s_adopted <= s_max,
        )

    # ------------------------------------------------------------------ #
    #  Affichage
    # ------------------------------------------------------------------ #

    def print_summary(self, n_legs: int = 2, cot_theta: float | None = None):
        """Imprime un résumé formaté du dimensionnement."""
        r = self.design_summary(n_legs, cot_theta)
        n_bars_label = f"{n_legs}×ø{self.phi_t*1e3:.0f}"
        print("\n=== Dimensionnement Effort Tranchant (EC2 §6.2.3) ===")
        print(f"  cot θ         = {r['cot_theta']:.2f}  (θ = {r['theta_deg']:.1f}°)")
        print(f"  V_Ed,max      = {r['v_ed_max']/1e3:.1f} kN")
        print(f"  V_Rd,max      = {r['v_rd_max']/1e3:.1f} kN  "
              f"{'✓' if r['ok_strut'] else '✗  ← augmenter bw ou f_ck'}")
        print(f"  Asw/st requis = {r['asw_s_max']*1e4:.3f} cm²/m  (max le long de la poutre)")
        print(f"  s requis      = {r['s_required_min']*1e2:.1f} cm  (espacement minimal)")
        print(f"  s_max (EC2)   = {r['s_max']*1e2:.1f} cm")
        print(f"  s adopté      = {r['s_adopted']*1e2:.0f} cm  "
              f"{'✓' if r['ok_spacing'] else '✗'}")
        print(f"  Cadres {n_bars_label} mm  "
              f"→ Asw/st fourni = {r['asw_s_provided']*1e4:.3f} cm²/m")