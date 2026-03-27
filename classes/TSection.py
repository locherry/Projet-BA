class TSection:
    def __init__(self, h_tot: float, b_eff: float, h_f: float, b_w: float) -> None:
        self.h_tot = h_tot
        self.b_eff = b_eff
        self.h_f = h_f
        self.b_w = b_w

    @property
    def h_w(self) -> float:
        return self.h_tot - self.h_f

    def areas(self) -> tuple[float, float]:
        A_f = self.b_eff * self.h_f
        A_w = self.b_w * self.h_w
        return A_f, A_w

    def centroid(self) -> float:
        A_f, A_w = self.areas()

        y_f = self.h_w + self.h_f / 2
        y_w = self.h_w / 2

        y_G = (A_f * y_f + A_w * y_w) / (A_f + A_w)
        return y_G

    def inertia_z(self) -> float:
        A_f, A_w = self.areas()
        y_G = self.centroid()

        y_f = self.h_w + self.h_f / 2
        y_w = self.h_w / 2

        I_f = (self.b_eff * self.h_f**3)/12 + A_f * (y_f - y_G)**2
        I_w = (self.b_w * self.h_w**3)/12 + A_w * (y_w - y_G)**2

        return I_f + I_w