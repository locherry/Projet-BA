from __future__ import annotations

from typing import Dict

from .Materials import Concrete, Steel
from .TSection import TSection
from .ReinforcementLayout import ReinforcementLayout
import numpy as np

from scipy.optimize import fsolve 


class FlexuralDesign:
    def __init__(
        self, section: TSection, concrete: Concrete, steel: Steel, d: float
    ) -> None:
        self.section = section
        self.concrete = concrete
        self.steel = steel
        self.d = d
        self.layout: ReinforcementLayout | None = None

    def neutral_axis(self, M_Ed: float) -> Dict[str, float]:
        b_eff = self.section.b_eff
        h_f = self.section.h_f
        f_cd = self.concrete.f_cd

        M_ct = b_eff * h_f * f_cd * (self.d - h_f / 2)

        if M_Ed < M_ct:
            mu = M_Ed / (b_eff * self.d**2 * f_cd)
            alpha = (1 - np.sqrt(1 - 2 * mu)) / 0.8
            x_na = alpha * self.d
            return {"x_na": x_na, "alpha": alpha, "mu": mu}
        else:
            raise NotImplementedError("Neutral axis in web not implemented")

    def reinforcement_area(self, M_Ed: float) -> float:
        res = self.neutral_axis(M_Ed)
        alpha = res["alpha"]

        epsilon_c = 3.5e-3
        epsilon_s = epsilon_c * (1 - alpha) / alpha
        print(epsilon_s)
        sigma_s = self.steel.stress(epsilon_s)

        F_c = 0.8 * alpha * self.section.b_eff * self.d * self.concrete.f_cd
        A_s = F_c / sigma_s

        return A_s

    def effective_concrete_modulus(self, phi: float = 2.5) -> float:
        """
        Effective modulus Ec,eff according to EC2 (creep coefficient phi).
        """
        return self.concrete.E_cm / (1.0 + phi)

    def steel_concrete_equivalence_ratio(self, phi: float = 2.5) -> float:
        """Equivalent n = E_s / E_c,eff."""
        return self.steel.Es / self.effective_concrete_modulus(phi)

    def phi_effective(self, phi_inf_t0: float, M_Eqp: float, M_Ed: float) -> float:
        """Effective creep coefficient for SLS: phi_ef = phi(inf,t0) * M_Eqp / M_Ed."""
        if M_Ed == 0:
            raise ValueError("M_Ed must be non-zero to compute phi_ef")
        return phi_inf_t0 * (M_Eqp / M_Ed)

    def check_stress_limitation(
        self,
        M_Ed: float,
        M_Eqp: float,
        # phi_inf_t0: float = 2.5,
        t0: float = 10,
    ) -> dict:
        if self.layout == None:
            raise(ValueError("Layout must be set"))
        """SLS stress check (compression and tension) based on EC2 rules and fixed n."""
        # Calculation phi_inf_t0
        a1 = (35 / (self.concrete.f_cm * 1e-6)) ** 0.7
        a2 = (35 / (self.concrete.f_cm * 1e-6)) ** 0.2
        h0 = (
            (self.section.b_w * self.section.h_tot)
            / (self.section.b_w + self.section.h_tot)
            * 1e3
        )
        RH = 50 / 100
        phiRH = (1 + (1 - RH) / (0.1 * h0 ** (1 / 3)) * a1) * a2
        betafcm = 16.8 / np.sqrt(self.concrete.f_cm * 1e-6)
        betat0 = 1 / (0.1 + t0**0.2)
        phi_inf_t0 = phiRH * betafcm * betat0

        # Current effective phi from SLS combination
        phi_ef = self.phi_effective(phi_inf_t0, M_Eqp, M_Ed)
        E_c_eff = self.concrete.E_cm / (1.0 + phi_ef)

        n = self.steel.Es / E_c_eff

        # Stress
        d = self.layout.compute_d_reel()
        A_s = self.reinforcement_area(M_Ed)
        def func(x):
            return self.section.b_w*x**2/2+(self.section.b_eff - self.section.b_w)*self.section.h_f*(x-self.section.h_f/2)-n*A_s*(d-x)
        x = fsolve(func, 0)[0]
    
        if x >= self.section.h_f :
            I_hr = self.section.b_w*x**3/3 + (self.section.b_eff - self.section.b_w)*self.section.h_f**3/12+ (self.section.b_eff - self.section.b_w)*self.section.h_f*(x-self.section.h_f/2)**2 + n*A_s*(d-x)**2
        else :
            raise(ValueError("ta pas la bonne taille de section"))
        sigma_s = n*M_Ed*(d-x)/I_hr
        sigma_c = M_Ed*x/I_hr

        limit_compression = 0.45 * self.concrete.f_ck
        limit_tension = self.concrete.f_ctm

        return {
            "M_Ed": M_Ed,
            "M_Eqp": M_Eqp,
            "phi_inf_t0": phi_inf_t0,
            "phi_ef": phi_ef,
            "E_cm": self.concrete.E_cm,
            "E_c_eff": E_c_eff,
            "n": self.steel.Es / E_c_eff,
            "ok_n": abs(self.steel.Es / E_c_eff - n) < 1e-6,
            "limit_t": limit_tension,
            "t0": t0,
            "sigma_s" : sigma_s,
            "sigma_c" : sigma_c,
            "limites de contrainte respectée acier, béton" : (sigma_s<=limit_tension, sigma_c<=limit_compression)
        }

    def crack_control(self, M_Ed, M_Eqp, A_s):
        # ------------------------------------------------------------------ #
        #  k(x) – non-linear stress distribution factor                      #
        # ------------------------------------------------------------------ #
        def k(x):
            """repartition non linéaire des contraintes de traction fct sur h (section
            rectangulaire ou section en T avec l'âme tendu) ou beff (section en T avec la
            membrane tendue)
            REF : slideshow 7, slide 20
            """
            if x < self.section.h_f:
                var = self.section.b_eff
            else:
                var = self.section.h_tot

            if var < 0.300:
                return 1
            elif var >= 0.800:
                return 0.65
            else:
                return 1.21 - 0.7 * var
        
        # ------------------------------------------------------------------ #
        #  1. Minimum reinforcement area  (slide 20)                          #
        # ------------------------------------------------------------------ #
        k_c = 0.4 #Flexion simple 
        

        # Concrete Area under tension before first crack
        x_na = self.neutral_axis(M_Ed)["x_na"]
        A_ct = self.section.b_w * (self.section.h_tot - x_na)
        
        # Maximum effort in tensile zone
        f_ct_eff = self.concrete.f_ctm  # f_ct,eff = f_ct,m  (slide 20)

        As_min = k_c * k(x_na) * f_ct_eff * A_ct / self.steel.f_yk
        
        print(As_min)
