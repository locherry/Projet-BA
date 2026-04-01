import numpy as np

class DistributedLoad:
    def __init__(self, q: float) -> None:
        self.q = q  # N/m


class Beam:
    def __init__(self, L: float, load: DistributedLoad) -> None:
        self.L = L
        self.load = load

    def reactions(self) -> tuple[float, float]:
        q, L = self.load.q, self.L
        Ra = q * L / 2
        Rb = q * L / 2
        return Ra, Rb

    def shear(self, x: np.ndarray) -> np.ndarray:
        Ra, _ = self.reactions()
        return Ra - self.load.q * x

    def moment(self, x: np.ndarray) -> np.ndarray:
        Ra, _ = self.reactions()
        return Ra * x - self.load.q * x**2 / 2

    def max_moment(self, q: float | None = None) -> float:
        if q is None:
            q = self.load.q
        L = self.L
        return (q * L**2) / 8