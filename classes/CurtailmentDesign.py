"""
CurtailmentDesign.py
====================
Épure d'arrêt des barres — Bar Curtailment Design  (EC2 §9.2.1.3)

Methodology
-----------
1. Plot the ULS bending moment envelope M_Ed(x)
2. Shift M_Ed by a_l = z·cotθ/2 toward the supports (unfavourable shift)
3. Compute M_Rd for sections with 1, 2, … N layers of longitudinal reinforcement
4. Find cut-off points where M_Rd(k layers) = M_Ed,shifted(x)
5. Verify anchorage length l_bd beyond each theoretical cut-off point (EC2 §8.4)

Usage (drop-in after ShearDesign in main.py)
--------------------------------------------
    from classes.CurtailmentDesign import CurtailmentDesign

    curtailment = CurtailmentDesign(
        beam=beam_uls,
        section=section,
        design=design,      # FlexuralDesign instance (already has layout attached)
        concrete=concrete,
        steel=steel,
        cot_theta=2.5,      # from ShearDesign
        x=x,                # same abscissa array used for shear
    )
    curtailment.print_summary()
    curtailment.plot()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
#  Data containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class LayerInfo:
    """One horizontal row of longitudinal reinforcement."""
    label: str          # e.g. "4HA32"
    n_bars: int
    phi: float          # bar diameter [m]
    As: float           # total area for this layer [m²]


@dataclass
class CutoffPoint:
    """Theoretical cut-off for one layer (left and right side)."""
    layer_index: int    # 0-based index of the layer being dropped
    x_left: float       # cut-off abscissa on left side [m]
    x_right: float      # cut-off abscissa on right side [m]
    lbd: float          # required anchorage length [m]


# ─────────────────────────────────────────────────────────────────────────────
#  CurtailmentDesign
# ─────────────────────────────────────────────────────────────────────────────

class CurtailmentDesign:
    """
    Bar curtailment (épure d'arrêt) for a simply-supported beam.

    Parameters
    ----------
    beam      : Beam instance (provides .moment(x) and .shear(x))
    section   : TSection instance
    design    : FlexuralDesign instance with layout already attached
    concrete  : Concrete instance  (needs .f_ck, .f_cd = f_ck/gamma_c)
    steel     : Steel instance     (needs .f_yk, .f_yd = f_yk/gamma_s)
    cot_theta : bielle angle cot θ from ShearDesign (default 2.5)
    x         : abscissa array [m] — reuse the one from ShearDesign
    n_points  : resolution for intersection search (ignored if x is provided)
    """

    def __init__(
        self,
        beam,
        section,
        design,
        concrete,
        steel,
        cot_theta: float = 2.5,
        x: np.ndarray | None = None,
        n_points: int = 2001,
    ):
        self.beam      = beam
        self.section   = section
        self.design    = design
        self.concrete  = concrete
        self.steel     = steel
        self.cot_theta = float(cot_theta)

        self.L   = beam.L
        self.d   = design.d          # effective depth [m]  (d_reel)
        self.z   = 0.9 * self.d      # lever arm [m]
        self.f_cd = concrete.f_ck / 1.5
        self.f_yd = steel.f_yk / 1.15

        self.x = (np.asarray(x, dtype=float) if x is not None
                  else np.linspace(0.0, self.L, n_points))

        # a_l shift: EC2 §9.2.1.3(2) — members with shear reinforcement
        #   a_l = z·cotθ / 2
        self.al = self.z * self.cot_theta / 2.0

        # Build layer list from the ReinforcementLayout attached to design
        self.layers: List[LayerInfo] = self._extract_layers()

        # Compute M_Ed and shifted M_Ed
        self.m_ed         = self.beam.moment(self.x)
        self.m_ed_shifted = self._shift_m_ed()

        # M_Rd for cumulative layer subsets  (index k → layers 0…k active)
        self.m_rd_cumul: List[float] = self._compute_m_rd_cumul()

        # Cut-off points
        self.cutoffs: List[CutoffPoint] = self._find_cutoffs()

    # ------------------------------------------------------------------ #
    #  Layer extraction from ReinforcementLayout
    # ------------------------------------------------------------------ #

    def _extract_layers(self) -> List[LayerInfo]:
        """
        Convert ReinforcementLayout rows into LayerInfo objects.
        Each row in the layout becomes one layer for curtailment purposes.
        """
        layout = self.design.layout
        layers = []
        for row in layout.rows:
            n_bars = row["n_groups"] * row["bars_per_group"]
            phi    = row["phi"]
            As     = n_bars * np.pi * (phi / 2) ** 2
            phi_mm = int(round(phi * 1e3))
            label  = f"{n_bars}HA{phi_mm}"
            layers.append(LayerInfo(label=label, n_bars=n_bars, phi=phi, As=As))
        return layers

    # ------------------------------------------------------------------ #
    #  Shifted M_Ed
    # ------------------------------------------------------------------ #

    def _shift_m_ed(self) -> np.ndarray:
        """
        Unfavourable shift of M_Ed by a_l toward the nearest support.

        For each abscissa x, the shifted demand equals M_Ed evaluated at
        a point moved a_l toward the nearer support:
            x_eff = x − a_l  (left half)
            x_eff = x + a_l  (right half)
        then clamped to [0, L].
        """
        # x_eff = np.where(
        #     self.x <= self.L / 2,
        #     np.maximum(0.0, self.x - self.al),
        #     np.minimum(self.L,  self.x + self.al),
        # )
        x_eff = np.where(
            self.x <= self.L / 2,
            self.x + self.al,
            self.x - self.al,
        )
        return self.beam.moment(x_eff)

    # ------------------------------------------------------------------ #
    #  M_Rd for cumulative subsets of layers
    # ------------------------------------------------------------------ #

    def _compute_m_rd_cumul(self) -> List[float]:
        """
        M_Rd for sections carrying 1, 2, … N cumulative layers [N·m].

        Uses the EC2 rectangular stress block (simplified method):
        - If N.A. is within the flange  → pure rectangular b_eff section
        - If N.A. is in the web         → composite T-section
        """
        b_eff = self.section.b_eff
        h_f   = self.section.h_f
        b_w   = self.section.b_w
        d     = self.d
        f_cd  = self.f_cd
        f_yd  = self.f_yd

        m_rd_list = []
        As_cumul  = 0.0

        for layer in self.layers:
            As_cumul += layer.As
            T = As_cumul * f_yd

            # Neutral axis depth assuming full flange width
            x_na = T / (0.8 * b_eff * f_cd)

            if x_na <= h_f:
                # N.A. in flange
                Mrd = T * (d - 0.4 * x_na)
            else:
                # N.A. in web
                C_fl  = 0.8 * (b_eff - b_w) * h_f * f_cd
                x_web = (T - C_fl) / (0.8 * b_w * f_cd)
                Mrd   = C_fl * (d - h_f / 2) + 0.8 * b_w * x_web * f_cd * (d - 0.4 * x_web)

            m_rd_list.append(float(Mrd))

        return m_rd_list

    # ------------------------------------------------------------------ #
    #  Anchorage length l_bd  (EC2 §8.4)
    # ------------------------------------------------------------------ #

    def _l_bd(self, phi: float) -> float:
        """
        Design anchorage length l_bd [m] for a bar of diameter phi [m].

        EC2 §8.4.2 :  f_bd = 2.25 · η₁ · η₂ · f_ctd
          η₁ = 1.0  (good bond conditions — bottom bars cast horizontally)
          η₂ = 1.0  (φ ≤ 32 mm)
        f_ctd = f_ctk,0.05 / γ_c = 0.7·0.3·f_ck^(2/3) / 1.5

        EC2 §8.4.3 :  l_b,rqd = (φ/4) · (σ_sd / f_bd)
        EC2 §8.4.4 :  l_bd = α₁·α₂·α₃·α₄·α₅ · l_b,rqd   (all α = 1 → conservative)
                      l_bd ≥ max(0.3·l_b,rqd ; 10·φ ; 100 mm)
        """
        f_ck_MPa = self.concrete.f_ck / 1e6
        f_ctd    = 0.7 * 0.3 * f_ck_MPa ** (2 / 3) / 1.5   # [MPa]
        f_bd     = 2.25 * 1.0 * 1.0 * f_ctd    
        
        # σ_sd at the cut-off is taken conservatively as f_yd (bar fully stressed).
        # EC2 §8.4.3 allows using the actual stress at the section, which would
        # reduce l_bd where bars are not fully yielding — conservative as-is.
        # [MPa]
        sigma_sd = self.f_yd / 1e6                           # [MPa]
        phi_mm   = phi * 1e3                                 # [mm]

        l_b_rqd = (phi_mm / 4) * (sigma_sd / f_bd)          # [mm]
        l_bd_mm = max(0.3 * l_b_rqd, 10 * phi_mm, 100.0)    # [mm]
        return l_bd_mm * 1e-3                                 # [m]

    # ------------------------------------------------------------------ #
    #  Intersection search
    # ------------------------------------------------------------------ #

    def _find_cutoffs(self) -> List[CutoffPoint]:
        """
        For each layer k (0-based), find the abscissae where
        M_Ed_shifted(x) crosses M_Rd[k] on the left and right halves.
        Layer k is theoretically unnecessary outside [x_left, x_right].
        """
        cutoffs = []
        x  = self.x
        ms = self.m_ed_shifted
        L  = self.L

        for k, mrd in enumerate(self.m_rd_cumul):
            lbd = self._l_bd(self.layers[k].phi)

            # --- left intersection (x < L/2, M increasing) ---
            x_left = self._interpolate_crossing(x, ms, mrd, side="left")

            # --- right intersection (x > L/2, M decreasing) ---
            x_right = self._interpolate_crossing(x, ms, mrd, side="right")

            cutoffs.append(CutoffPoint(
                layer_index=k,
                x_left=x_left,
                x_right=x_right,
                lbd=lbd,
            ))

        return cutoffs

    @staticmethod
    def _interpolate_crossing(
        x: np.ndarray,
        m: np.ndarray,
        m_target: float,
        side: str,
    ) -> float:
        """Linear interpolation of the abscissa where m(x) = m_target."""
        L = x[-1]
        if side == "left":
            xs, ms = x[x <= L / 2], m[x <= L / 2]
            # ascending crossing
            for i in range(1, len(xs)):
                if ms[i - 1] < m_target <= ms[i]:
                    t = (m_target - ms[i - 1]) / (ms[i] - ms[i - 1])
                    return float(xs[i - 1] + t * (xs[i] - xs[i - 1]))
            return 0.0   # M_Rd never reached → bar needed from support
        else:
            xs, ms = x[x >= L / 2], m[x >= L / 2]
            # descending crossing
            for i in range(1, len(xs)):
                if ms[i - 1] >= m_target > ms[i]:
                    t = (m_target - ms[i - 1]) / (ms[i] - ms[i - 1])
                    return float(xs[i - 1] + t * (xs[i] - xs[i - 1]))
            return L     # M_Rd never reached → bar needed to support

    # ------------------------------------------------------------------ #
    #  Summary
    # ------------------------------------------------------------------ #

    def summary(self) -> List[dict]:
        """
        Return a list of dicts, one per layer, with curtailment results.
        """
        results = []
        As_cumul = 0.0
        for k, (layer, c) in enumerate(zip(self.layers, self.cutoffs)):
            As_cumul += layer.As
            results.append(dict(
                layer         = layer.label,
                As_layer      = layer.As,
                As_cumul      = As_cumul,
                M_Rd          = self.m_rd_cumul[k],
                x_cut_left    = c.x_left,
                x_cut_right   = c.x_right,
                bar_start     = max(0.0, c.x_left  - c.lbd),
                bar_end       = min(self.L, c.x_right + c.lbd),
                lbd           = c.lbd,
                ok_anchorage  = c.lbd <= c.x_left,   # enough room from support
            ))
        return results

    def print_summary(self):
        """Print a formatted curtailment summary to stdout."""
        print("\n" + "=" * 65)
        print("  ÉPURE D'ARRÊT DES BARRES  (EC2 §9.2.1.3 & §8.4)")
        print("=" * 65)
        print(f"  d       = {self.d * 1e2:.2f} cm")
        print(f"  z       = {self.z * 1e2:.2f} cm")
        print(f"  cotθ    = {self.cot_theta:.1f}")
        print(f"  a_l     = z·cotθ/2 = {self.al * 1e2:.2f} cm")
        print("-" * 65)

        for r in self.summary():
            print(f"\n  Layer : {r['layer']}")
            print(f"    As (layer / cumul)  = {r['As_layer']*1e4:.2f} / {r['As_cumul']*1e4:.2f} cm²")
            print(f"    M_Rd               = {r['M_Rd']/1e6:.3f} MN·m")
            print(f"    x_cut (L / R)      = {r['x_cut_left']:.3f} m  /  {r['x_cut_right']:.3f} m")
            print(f"    l_bd               = {r['lbd']*1e2:.1f} cm")
            print(f"    Bar extent         = [{r['bar_start']:.3f} ; {r['bar_end']:.3f}] m")
            ok = "✓" if r["ok_anchorage"] else "✗  (check end anchorage at support)"
            print(f"    Anchorage check    = {ok}")

        print("=" * 65)