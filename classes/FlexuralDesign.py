from __future__ import annotations

from typing import Dict

from .Materials import Concrete, Steel
from .TSection import TSection
from .ReinforcementLayout import ReinforcementLayout
import numpy as np


class FlexuralDesign:
    def __init__(self, section: TSection, concrete: Concrete, steel: Steel, d: float) -> None:
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
        phi_inf_t0: float = 2.5,
        t0: float = 10.0,
    ) -> dict:
        """SLS stress check (compression and tension) based on EC2 rules and fixed n."""
        # Current effective phi from SLS combination
        phi_ef = self.phi_effective(phi_inf_t0, M_Eqp, M_Ed)
        E_eff = self.concrete.E_cm / (1.0 + phi_ef)

        limit_compression = 0.45 * self.concrete.f_ck
        limit_tension = self.concrete.f_ctm
        
        n = self.steel.Es / E_eff

        return {
            "M_Ed": M_Ed,
            "M_Eqp": M_Eqp,
            "phi_inf_t0": phi_inf_t0,
            "phi_ef": phi_ef,
            "E_cm": self.concrete.E_cm,
            "E_eff": E_eff,
            "n": self.steel.Es / E_eff,
            "ok_n": abs(self.steel.Es / E_eff - n) < 1e-6,
            "limit_t": limit_tension,
            "t0": t0,
        }
