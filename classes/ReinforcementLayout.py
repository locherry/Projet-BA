from __future__ import annotations
from typing import List, Tuple
from .TSection import TSection
import numpy as np


class ReinforcementLayout:
    """
    Represents the longitudinal reinforcement layout of a T-section.
    Computes effective depth (d_reel) and provides bar positions for plotting.
    Supports vertically stacked groups.
    """

    def __init__(self, section: TSection, phi_t: float = 0.010):
        """
        Parameters
        ----------
        section : TSection
        phi_t   : float
            Stirrup (link) diameter in metres. Must be provided explicitly —
            e.g. 0.010 for ø10 mm.  Default 10 mm is a safe starting value;
            always override it with the actual design choice.
        """
        self.section = section
        self.rows: List[dict] = []   # ordered bottom → top
        self._phi_t = phi_t          # stirrup diameter is a design input, not derived

    # ---------------------------------------------------------------------- #
    #  Row management                                                          #
    # ---------------------------------------------------------------------- #

    def add_row(
        self,
        n_groups: int,
        bars_per_group: int,
        diameter: float,
        grouped: bool = False,
    ):
        """
        Add a horizontal row of bars (or bar groups).

        Parameters
        ----------
        n_groups       : number of groups (or single bars) across the width
        bars_per_group : bars bundled vertically in each group
        diameter       : individual bar diameter [m]
        grouped        : True → bars are stacked vertically (bundle)
        """
        self.rows.append({
            "n_groups": n_groups,
            "bars_per_group": bars_per_group,
            "phi": diameter,
            "grouped": grouped,
        })

    # ---------------------------------------------------------------------- #
    #  Helper diameters                                                        #
    # ---------------------------------------------------------------------- #

    def equivalent_diameter(self, row: dict) -> float:
        """Equivalent diameter of a group: φ_n = φ · √n  (EC2 §8.9.1)."""
        if row["grouped"]:
            return row["phi"] * np.sqrt(row["bars_per_group"])
        return row["phi"]

    def phi_eq(self) -> float:
        """
        Global equivalent bar diameter for crack-spacing formula.
        φ_eq = Σ(nᵢ φᵢ²) / Σ(nᵢ φᵢ)   (EC2 §7.3.4(3))
        """
        sum_n_phi2 = 0.0
        sum_n_phi = 0.0
        for row in self.rows:
            n_i = row["n_groups"] * row["bars_per_group"]
            phi_i = row["phi"]
            sum_n_phi2 += n_i * phi_i ** 2
            sum_n_phi += n_i * phi_i
        if sum_n_phi == 0:
            raise ValueError("No reinforcement rows added yet")
        return sum_n_phi2 / sum_n_phi

    # ---------------------------------------------------------------------- #
    #  Cover                                                                   #
    # ---------------------------------------------------------------------- #

    @property
    def c_nom(self) -> float:
        """
        Nominal concrete cover  c_nom = c_min + Δc_dev  (EC2 §4.4.1).

        c_min = max(c_min,b ; c_min,dur ; 10 mm)

        c_min,b  = φ_eq,max           (EC2 §4.4.1.2(1) — single bar)
                 = φ_n  for bundles   (EC2 §4.4.1.2(3))
        c_min,dur is taken from the exposure-class parameter supplied at
                  construction (default XC1 → 15 mm for reinforced concrete,
                  EC2 Table 4.4N).
        Δc_dev   = 10 mm              (EC2 §4.4.1.3, recommended value)
        """
        if not self.rows:
            raise ValueError("No reinforcement rows added yet")
        # largest equivalent diameter governs bond cover
        phi_eq_max = max(self.equivalent_diameter(r) for r in self.rows)
        c_min_b = phi_eq_max                # EC2 §4.4.1.2 — no extra 5 mm
        c_min_dur = self._c_min_dur         # from exposure class
        delta_c_dev = 10e-3                 # EC2 §4.4.1.3
        c_min = max(c_min_b, c_min_dur, 10e-3)
        return c_min + delta_c_dev

    # Exposure-class driven durability cover — override via subclass or setter
    # Default: XC1, structural class S4 → c_min,dur = 15 mm (EC2 Table 4.4N)
    _c_min_dur: float = 15e-3

    def set_exposure(self, c_min_dur: float):
        """Set c_min,dur [m] according to EC2 Table 4.4N for the actual exposure class."""
        self._c_min_dur = c_min_dur

    # ---------------------------------------------------------------------- #
    #  Stirrup diameter (design input)                                         #
    # ---------------------------------------------------------------------- #

    @property
    def phi_t(self) -> float:
        """Stirrup diameter [m] — supplied at construction, not derived."""
        return self._phi_t

    # ---------------------------------------------------------------------- #
    #  Effective depth  d_reel                                                 #
    # ---------------------------------------------------------------------- #

    def compute_d_reel(self) -> float:
        """
        Compute the effective depth as the area-weighted centroid of ALL bars
        defined via add_row(), measured from the bottom fibre.

        Uses the actual bar positions returned by get_bar_positions() so that
        d_reel is always consistent with the real layout — no hardcoded values.
        """
        if not self.rows:
            raise ValueError("No reinforcement rows added yet")

        total_area = 0.0
        weighted_y = 0.0

        for x, y, r in self.get_bar_positions():
            A_bar = np.pi * r ** 2          # area of one individual bar
            total_area += A_bar
            weighted_y += A_bar * y

        if total_area == 0:
            raise ValueError("Total reinforcement area is zero")

        d_reel = self.section.h_tot - weighted_y / total_area
        return d_reel

    # ---------------------------------------------------------------------- #
    #  Bar positions for plotting                                              #
    # ---------------------------------------------------------------------- #

    def get_bar_positions(self) -> List[Tuple[float, float, float]]:
        """
        Return a list of (x, y, radius) for every individual bar,
        where y is measured from the **bottom** fibre.

        Rows are ordered bottom → top; within each row, groups are
        spread evenly across the web width and stacked vertically when grouped.
        """
        positions: List[Tuple[float, float, float]] = []
        d_g = 12e-3          # maximum aggregate size (assumed 12 mm)
        y_current = 0.0      # centroid of current stack from bottom fibre

        for i, row in enumerate(self.rows):
            phi = row["phi"]
            n_groups = row["n_groups"]
            bars_per_group = row["bars_per_group"] if row["grouped"] else 1
            stack_height = phi * bars_per_group

            # ---- vertical position of stack centroid ----
            if i == 0:
                # First (bottom) row: cover + stirrup + half stack
                y_current = self.c_nom + self.phi_t + stack_height / 2
            else:
                prev_row = self.rows[i - 1]
                phi_prev = prev_row["phi"]
                bars_prev = prev_row["bars_per_group"] if prev_row["grouped"] else 1
                stack_prev = phi_prev * bars_prev
                # Clear spacing ≥ max(φ, φ_prev, d_g + 5 mm, 20 mm)  EC2 §8.2
                e_v = max(phi, phi_prev, d_g + 5e-3, 20e-3)
                y_current += stack_prev / 2 + e_v + stack_height / 2

            # ---- horizontal positions (centred on web axis) ----
            b_w = self.section.b_w
            cover = self.c_nom + self.phi_t
            usable_width = b_w - 2 * cover
            x_positions = (
                np.linspace(-usable_width / 2, usable_width / 2, n_groups)
                if n_groups > 1
                else [0.0]
            )

            # ---- place individual bars within each group stack ----
            for xg in x_positions:
                for k in range(bars_per_group):
                    y_bar = y_current - stack_height / 2 + phi / 2 + k * phi
                    positions.append((float(xg), float(y_bar), phi / 2))

        return positions