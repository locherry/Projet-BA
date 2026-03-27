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
        self.layout: ReinforcementLayout | None = None  # optional typing

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
