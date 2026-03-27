class Concrete:
    def __init__(self, f_ck: float, gamma_c: float) -> None:
        self.f_ck = f_ck
        self.gamma_c = gamma_c

    @property
    def f_cd(self) -> float:
        return self.f_ck / self.gamma_c


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