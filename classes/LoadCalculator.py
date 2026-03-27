class LoadCalculator:
    def __init__(self) -> None:
        # Dead loads (G) in N/m
        self.G_components = {
            'charpente': 1e3,
            'bac': 4e3,
            'rev1': 2.4e3,
            'voile': 0.6e3,
            'dalle': 7.5e3,
            'retombee': 13.75e3,
            'rev2': 5.6e3,
            'scenographie': 16e3,
        }
        # Live load in N/m
        self.Q_components = {
            'bureaux': 5e3,
            'salle': 10e3
        }
    
    @property
    def G_total(self) -> float:
        return sum(self.G_components.values())
    
    @property
    def Q_total(self) -> float:
        return sum(self.Q_components.values())
    
    @property
    def q_ELS(self) -> float:
        return self.G_total + self.Q_total
    
    @property
    def q_ELU(self) -> float:
        return 1.35 * self.G_total + 1.5 * self.Q_total