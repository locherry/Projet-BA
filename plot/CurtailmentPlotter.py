import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from typing import List, Tuple
import numpy as np

class CurtailmentPlotter:
    def __init__(self, curtailment, design):
        self.curtailment = curtailment
        self.design = design # reinforcement design
    def plot(
        self,
        figsize: Tuple[float, float] = (10, 5),
        save_path="./plot/curtailment.svg",
    ):
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
        M_ed    = self.curtailment.m_ed / 1e6
        M_shift = np.maximum(self.curtailment.m_ed_shifted / 1e6, M_ed)

        ax1.plot(x, M_ed,    color="tab:blue",   label="M_Ed")
        ax1.plot(x, M_shift, color="tab:blue", linestyle="--",
                label=f"M_Ed décalé (a_l={self.curtailment.al:.2f} m)")

        # ── M_Rd horizontal levels + cut-off verticals ──
        colors = ["tab:orange", "tab:green", "tab:red", "tab:purple"]
        for k, (c, mrd) in enumerate(
            zip(self.curtailment.cutoffs, self.curtailment.m_rd_cumul)
        ):
            col = colors[k % len(colors)]
            ax1.axhline(mrd / 1e6, linestyle=":", color=col,
                        label=f"M_Rd,{k+1} = {mrd/1e6:.3f} MN·m")
            if c.x_left > 0:
                ax1.axvline(c.x_left,  linestyle="--", color="tab:grey", linewidth=0.8)
                ax1.axvline(c.x_right, linestyle="--", color="tab:grey", linewidth=0.8)

        # ── Resistance envelope M_Rd(x) ──
        # Build piecewise:
        #   [0, x_cut_left[1] ]           → M_Rd[0]  (bottom layer alone, full length)
        #   [x_cut_left[1], x_cut_left[1]+lbd[1]] → ramp from M_Rd[0] to M_Rd[1]
        #   [x_cut_left[1]+lbd[1], L/2]   → M_Rd[1]  (all layers, midspan)
        #   symmetric on right half

        layers  = self.curtailment.layers
        cutoffs = self.curtailment.cutoffs
        m_rds   = self.curtailment.m_rd_cumul   # cumulative, index 0 = bottom layer only

        # Work on left half, then mirror
        x_env_l = [0.0]
        m_env_l = [m_rds[0] / 1e6]   # starts at M_Rd of bottom layer

        for k in range(1, len(layers)):
            c      = cutoffs[k]
            lbd    = c.lbd
            x_cut  = c.x_left          # theoretical cut-off on left
            m_low  = m_rds[k - 1] / 1e6
            m_high = m_rds[k]     / 1e6

            # flat at m_low until cut-off point
            x_env_l.append(x_cut)
            m_env_l.append(m_low)

            # ramp over l_bd
            x_env_l.append(x_cut + lbd)
            m_env_l.append(m_high)

        # extend flat to midspan
        x_env_l.append(L / 2)
        m_env_l.append(m_rds[-1] / 1e6)

        # mirror to right half
        x_env_r = [L - xi for xi in reversed(x_env_l)]
        m_env_r = list(reversed(m_env_l))

        x_env = x_env_l + x_env_r
        m_env = m_env_l + m_env_r

        ax1.plot(x_env, m_env, color="tab:red", linewidth=2.0,
                label="M_Rd(x) enveloppe")
        ax1.fill_between(x_env, m_env, alpha=0.08, color="tab:red")

        ax1.set_ylabel("M [MN·m]")
        ax1.set_title("Raccourcissement des barres longitudinales")
        ax1.grid(True)
        ax1.legend(fontsize=8)
        ax1.set_ylim(0, M_ed.max() * 1.2)

        # =========================
        # BOTTOM: REINFORCEMENT
        # =========================
        y = 0
        colors2 = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

        for k, (layer, cutoff) in enumerate(
            zip(self.curtailment.layers, self.curtailment.cutoffs)
        ):
            col = colors2[k % len(colors2)]
            xL  = cutoff.x_left
            xR  = cutoff.x_right
            lbd = cutoff.lbd

            # Main bar extent (cut-off to cut-off)
            ax2.plot([xL, xR], [y, y], linewidth=3, color=col, label=layer.label)

            # l_bd annotation (left side)
            if xL > 0:
                ax2.annotate(
                    "", xy=(xL, y + 0.3), xytext=(max(0, xL + lbd), y + 0.3),
                    arrowprops=dict(arrowstyle="<->", linewidth=0.7, color=col),
                )
                ax2.text((xL + max(0, xL + lbd)) / 2, y + 0.38,
                        f"l_bd={lbd*100:.0f}cm",
                        ha="center", fontsize=6.5, color=col)

            # Label
            ax2.text(xR + 0.01 * L, y, layer.label, va="center",
                    fontsize=8, color=col)

            y += 1.0

        ax2.set_xlabel("x [m]")
        ax2.set_yticks([])
        ax2.set_ylim(-0.5, y)
        ax2.grid(True, axis="x")
        ax2.set_xlim(0, L)

        plt.tight_layout()
        plt.savefig(save_path, dpi=160, bbox_inches="tight")
        print(f"Saved → {save_path}")