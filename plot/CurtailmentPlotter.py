import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from typing import List, Tuple


class CurtailmentPlotter:
    def __init__(self, curtailment, design):
        self.curtailment = curtailment
        self.design = design # reinforcement design

    def plot(
        self,
        figsize: Tuple[float, float] = (10, 5),
        save_path="./plot/curtailement.svg",
    ):
        """
        Top:
            - M_Ed
            - Shifted M_Ed
            - M_Rd levels
            - Cut-off points

        Bottom:
            - Reinforcement layout
            - Anchorage lengths (l_bd)

        Shared x-axis for consistency.
        """
        fig, (ax1, ax2) = plt.subplots(
            2, 1,
            figsize=figsize,
            sharex=True,
            gridspec_kw={"height_ratios": [2, 1]}
        )

        x = self.curtailment.x
        L = self.curtailment.L

        # =========================
        # TOP: MOMENT DIAGRAM
        # =========================
        M_ed = self.curtailment.m_ed / 1e6
        M_shift = self.curtailment.m_ed_shifted / 1e6

        ax1.plot(x, M_ed, label="M_Ed")
        ax1.plot(x, M_shift, "--", label=f"M_Ed shifted (a_l={self.curtailment.al:.2f} m)")

        for k, (c, mrd) in enumerate(
            zip(self.curtailment.cutoffs, self.curtailment.m_rd_cumul)
        ):
            mrd_val = mrd / 1e6

            ax1.axhline(mrd_val, linestyle=":", label=f"M_Rd,{k+1}")
            ax1.axvline(c.x_left, linestyle="--")
            ax1.axvline(c.x_right, linestyle="--")

        ax1.set_ylabel("M [MN·m]")
        ax1.set_title("Raccourcissement des barres longitudinales")
        ax1.grid(True)
        ax1.legend(fontsize=8)

        # =========================
        # BOTTOM: REINFORCEMENT
        # =========================
        y_step = 1.0
        y = 0

        for k, (layer, cutoff) in enumerate(
            zip(self.curtailment.layers, self.curtailment.cutoffs)
        ):
            xL = cutoff.x_left
            xR = cutoff.x_right

            # Full bar (anchorage included)
            ax2.plot([xL, xR], [y, y], linewidth=3)

            # Label
            ax2.text(
                xR + 0.01 * L,
                y,
                layer.label,
                va="center",
                fontsize=8,
            )

            y += y_step

        ax2.set_xlabel("x [m]")
        ax2.set_yticks([])
        ax2.set_ylim(-0.5, y)
        ax2.grid(True, axis="x")

        # =========================
        # FINALIZE
        # =========================
        ax2.set_xlim(0, L)

        plt.tight_layout()
        plt.savefig(save_path, dpi=160, bbox_inches="tight")
        print(f"Saved → {save_path}")