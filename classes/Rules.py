from dataclasses import dataclass

@dataclass(frozen=True)
class DurabilityRules:
    @staticmethod
    def structural_class(working_life_years: int) -> str:
        """
        Determine the structural class (S1 to S6) based on the working life in years.
        Ref : slideshow 5 slide 48
        """
        
        if working_life_years <= 25:
            return "S3"
        if working_life_years <= 50:
            return "S4"
        if working_life_years <= 100:
            return "S5"
        return "S6"

    @staticmethod
    def c_min_dur(exposure_class: str, working_life_years: int) -> float:
        """
        Compute the durability-driven minimum cover c_min.
        
        Ref : slideshow 5 slide 49
        """
        sc = DurabilityRules.structural_class(working_life_years)
        table = {
            "XC1": {"S1": 10e-3, "S2": 10e-3, "S3": 10e-3, "S4": 10e-3, "S5": 15e-3, "S6": 20e-3},
            "XC2": {"S1": 20e-3, "S2": 20e-3, "S3": 25e-3, "S4": 30e-3, "S5": 35e-3, "S6": 40e-3},
            "XC3": {"S1": 20e-3, "S2": 20e-3, "S3": 25e-3, "S4": 30e-3, "S5": 35e-3, "S6": 40e-3},
            "XC4": {"S1": 25e-3, "S2": 25e-3, "S3": 30e-3, "S4": 35e-3, "S5": 40e-3, "S6": 45e-3},
        }
        return table[exposure_class.upper()][sc]

    @staticmethod
    def w_max(exposure_class: str, load_combination: str = "quasi_permanent") -> float:
        """
        Return the EC2 recommended maximum crack width w_max in meters.
        REF : slideshow 6, slide 22

        Parameters
        ----------
        exposure_class : str
            Exposure class such as X0, XC1, XC2, XC3, XC4, XD1, XS1, ...
        load_combination : str
            "quasi_permanent" or "frequent" (for prestressed members / special cases)

        Returns
        -------
        float
            w_max in meters
        """
        ec = exposure_class.upper()
        if ec in {"X0", "XC1"}:
            return 0.4e-3
        if ec in {"XC2", "XC3", "XC4"}:
            return 0.3e-3
        if ec in {"XD1", "XD2", "XS1", "XS2", "XS3"}:
            return 0.3e-3 if load_combination == "quasi_permanent" else 0.2e-3
        return 0.3e-3