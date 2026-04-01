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

    def __init__(self, section: TSection):
        self.section = section
        self.rows: List[dict] = []  # bottom → top

    def add_row(self, n_groups: int, bars_per_group: int, diameter: float, grouped: bool = False):
        self.rows.append({
            "n_groups": n_groups,
            "bars_per_group": bars_per_group,
            "phi": diameter,
            "grouped": grouped
        })

    def equivalent_diameter(self, row: dict) -> float:
        if row["grouped"]:
            return row["phi"] * np.sqrt(row["bars_per_group"])
        return row["phi"]
    
    def phi_eq(self) -> float:
        """Equivalent bar diameter φ_eq = Σ(n_i φ_i²) / Σ(n_i φ_i)  (EC2 slide 29)."""
        sum_n_phi2 = 0.0
        sum_n_phi  = 0.0
        for row in self.rows:
            n_i   = row["n_groups"] * row["bars_per_group"]
            phi_i = row["phi"]
            sum_n_phi2 += n_i * phi_i**2
            sum_n_phi  += n_i * phi_i
        return sum_n_phi2 / sum_n_phi

    @property
    def c_nom(self) -> float:
        """Nominal concrete cover based on largest equivalent bar."""
        if not self.rows:
            raise ValueError("No reinforcement rows added yet")
        phi_eq_max = max(self.equivalent_diameter(r) for r in self.rows)
        c_min_b = phi_eq_max + 5e-3
        c_min_dur = 10e-3
        delta_c_dev = 10e-3
        return max(c_min_b, c_min_dur, 10e-3) + delta_c_dev

    @property
    def phi_t(self) -> float:
        """Stirrup diameter, 30% of the largest equivalent bar."""
        if not self.rows:
            raise ValueError("No reinforcement rows added yet")
        phi_eq_max = max(self.equivalent_diameter(r) for r in self.rows)
        return 0.3 * phi_eq_max

    def compute_d_reel(self) -> float:
        """
        Compute the effective depth of the reinforcement based on the fixed configuration:
        - Top row: 1 bar of 32 mm
        - Bottom row: 3 groups of 2 bars each (vertically stacked)
        Uses class properties c_nom and phi_t.
        """
        # Bar diameters
        phi_1 = 32e-3                  # Top single bar
        phi_2 = phi_1 * np.sqrt(2)     # Bottom groups of 2 bars

        # Steel areas
        A_s_1 = np.pi * (phi_1 / 2) ** 2
        A_s_2 = np.pi * (phi_2 / 2) ** 2

        # Use class properties
        c_nom = self.c_nom
        phi_t = self.phi_t

        # Vertical spacing
        d_g = 12e-3
        e_v_min = max(phi_2, d_g + 5e-3, 20e-3)
        e_v = e_v_min

        # Horizontal spacing check
        e_h = self.section.b_w - 2 * c_nom - phi_2 * 3 - 2 * phi_t
        if e_h < e_v_min:
            raise ValueError(f"e_h({e_h:.3f}) < e_v_min ({e_v:.3f})")

        # Effective depth per row
        d_2 = c_nom + phi_t + phi_2 / 2          # Bottom row centroid
        d_1 = c_nom + phi_t + phi_2 + e_v + phi_1 / 2  # Top row centroid

        # Weighted average for effective depth
        A_s_tot = A_s_1 * 1 + A_s_2 * 3
        d_reel = (d_2 * A_s_2 * 3 + d_1 * A_s_1 * 1) / A_s_tot

        return d_reel
    
    def get_bar_positions(self) -> List[Tuple[float, float, float]]:
        """
        Returns a list of (x, y, radius) for each individual bar.
        Groups are stacked vertically, and bars within a group are centered.
        """
        positions: List[Tuple[float, float, float]] = []
        d_g = 12e-3
        y_current = 0

        for i, row in enumerate(self.rows):
            phi = row["phi"]
            n_groups = row["n_groups"]
            bars_per_group = row["bars_per_group"] if row["grouped"] else 1
            stack_height = phi * bars_per_group

            # First row vertical position
            if i == 0:
                y_current = self.c_nom + self.phi_t + stack_height / 2
            else:
                prev_row = self.rows[i - 1]
                phi_prev = prev_row["phi"]
                bars_prev = prev_row["bars_per_group"] if prev_row["grouped"] else 1
                stack_prev = phi_prev * bars_prev
                e_v = max(phi, phi_prev, d_g + 5e-3, 20e-3)
                y_current += stack_prev / 2 + e_v + stack_height / 2

            # Horizontal positions
            b_w = self.section.b_w
            cover = self.c_nom + self.phi_t
            usable_width = b_w - 2 * cover
            x_positions = np.linspace(-usable_width / 2, usable_width / 2, n_groups) if n_groups > 1 else [0]

            # Place bars vertically in stack
            for xg in x_positions:
                for k in range(bars_per_group):
                    y_bar = y_current - stack_height / 2 + phi / 2 + k * phi
                    positions.append((xg, y_bar, phi / 2))

        return positions