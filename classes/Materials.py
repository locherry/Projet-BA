import numpy as np


class Concrete:
    def __init__(self, f_ck: float, gamma_c: float) -> None:
        self.f_ck = f_ck
        self.gamma_c = gamma_c

    @property
    def f_cd(self) -> float:
        return self.f_ck / self.gamma_c

    @property
    def f_cm(self) -> float:
        """
        mean concrete strength
        Ref :characteristic of concrete table : slideshow 3, slide 20
        """
        # mean cylinder strength, f_cm = f_ck + 8 MPa (Eurocode 2)
        return self.f_ck + 8e6

    @property
    def E_cm(self) -> float:
        """
        elastic modulus
        Ref :characteristic of concrete table : slideshow 3, slide 20
        """
        # EC2: E_cm = 22*(f_cm/10)^0.3 [GPa] -> convert to Pa
        return 22e9 * ((self.f_cm / 1e6) / 10) ** 0.3

    @property
    def f_ctm(self) -> float:
        """
        mean tensile strength
        Ref :characteristic of concrete table : slideshow 3, slide 20
        """
        if self.f_ck <= 50e6:
            # Mean tensile strength (EC2): 0.30 * f_ck^(2/3), f_ck in MPa
            # convert to Pa
            return 0.30 * (self.f_ck / 1e6) ** (2 / 3) * 1e6
        else:
            return 2.12 * np.log(1 + self.f_cm / 10)


class Steel:
    def __init__(self, f_yk: float, gamma_s: float, Es: float = 200e9) -> None:
        self.f_yk = f_yk
        self.gamma_s = gamma_s
        self.Es = Es

    @property
    def f_yd(self) -> float:
        return self.f_yk / self.gamma_s

    def stress(self, epsilon: float) -> float:
        epsilon_se = self.f_yd / self.Es

        if epsilon <= epsilon_se:
            return self.Es * epsilon
        else:
            return self.f_yd + 0.817e6 * (epsilon * 1000 - 2.17)
