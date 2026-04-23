import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class BeamPlotter:
    def __init__(self, beam):
        self.beam = beam

    # ── Support drawing ─────────────────────────────────────────────
    def _draw_pin(self, ax, xc, yc):
        h, w = 0.6, 0.5
        tri = plt.Polygon( # pyright: ignore[reportPrivateImportUsage]
            [[xc, yc], [xc - w/2, yc - h], [xc + w/2, yc - h]],
            closed=True, fc="lightgray", ec="k", lw=1.5, zorder=5
        )
        ax.add_patch(tri)

    def _draw_roller(self, ax, xc, yc):
        h = 0.5
        ax.add_patch(plt.Circle( # pyright: ignore[reportPrivateImportUsage]
            (xc, yc - h/2 - .15), h/2,
            fc="lightgray", ec="k", lw=1.5, zorder=5
        ))

    # ── Core computation (delegates to Beam) ────────────────────────
    def _compute_fields(self, n=1000):
        L = self.beam.L
        x = np.linspace(0, L, n)

        T = self.beam.shear(x)
        M = self.beam.moment(x)

        Ra, Rb = self.beam.reactions()

        T_max = Ra
        T_min = -Rb
        M_max = self.beam.max_moment()
        x_Mmax = L / 2

        return x, T, M, Ra, Rb, T_max, T_min, M_max, x_Mmax

    # ── Main plot ───────────────────────────────────────────────────
    def plot(self, save_path="./plot/beam_diagram.svg"):
        L = self.beam.L
        q = self.beam.load.q

        x, T, M, Ra, Rb, T_max, T_min, M_max, x_Mmax = self._compute_fields()

        # Convert all force/moment values from SI (N, N·m) to kN, kN·m
        q      = q      / 1e3
        T      = T      / 1e3
        M      = M      / 1e3
        Ra     = Ra     / 1e3
        Rb     = Rb     / 1e3
        T_max  = T_max  / 1e3
        T_min  = T_min  / 1e3
        M_max  = M_max  / 1e3

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
        fig.suptitle(
            "Poutre isostatique — TNP Villeurbanne  |  Niveau +7.55 m  |  "
            f"$L = {L}$ m  |  $q = {q}$ kN/m",
            fontsize=12, fontweight="bold"
        )
        fig.subplots_adjust(hspace=0.55, top=0.93, bottom=0.07)

        # ═══════════════════════════════════════════════════════════════
        # SUBPLOT 1 — Schéma statique
        # ═══════════════════════════════════════════════════════════════
        ax1.set_xlim(0, L)
        ax1.set_ylim(-2.6, 4.2)
        ax1.axis("off")
        ax1.set_title("Schéma statique", fontsize=11, fontweight="bold")

        beam_y, beam_h = 0.0, 0.30

        # Beam
        ax1.add_patch(mpatches.FancyBboxPatch(
            (0, beam_y - beam_h/2), L, beam_h,
            boxstyle="round,pad=0.04",
            fc="#D5D8DC", ec="k", lw=1.8, zorder=3
        ))

        # Distributed load
        arrow_top = 3.0
        arrow_bot = beam_y + beam_h / 2
        for xi in np.linspace(0.3, L - 0.3, 16):
            ax1.annotate("", xy=(xi, arrow_bot), xytext=(xi, arrow_top),
                         arrowprops=dict(arrowstyle="-|>", color="tab:red",
                                         lw=1.2, mutation_scale=8), zorder=4)
        ax1.plot([0, L], [arrow_top, arrow_top], lw=2.0, color="tab:red", zorder=4)
        ax1.text(L/2, arrow_top + 0.35, f"$q = {q}$ kN/m",
                 ha="center", va="bottom", fontsize=11,
                 color="tab:red", fontweight="bold")

        # Supports
        self._draw_pin(ax1, 0, beam_y - beam_h/2)
        self._draw_roller(ax1, L, beam_y - beam_h/2)

        # Reactions
        rl = 1.3
        for xc in [0, L]:
            ax1.annotate("", xy=(xc, beam_y - beam_h/2),
                         xytext=(xc, beam_y - beam_h/2 - rl),
                         arrowprops=dict(arrowstyle="-|>", color="tab:blue",
                                         lw=2.0, mutation_scale=14), zorder=6)

        ax1.text(0.22, beam_y - beam_h/2 - rl/2,
                 f"$R_E={Ra:.1f}$ kN", ha="left", va="center",
                 fontsize=10, color="tab:blue", fontweight="bold")

        ax1.text(L - 0.22, beam_y - beam_h/2 - rl/2,
                 f"$R_C={Rb:.1f}$ kN", ha="right", va="center",
                 fontsize=10, color="tab:blue", fontweight="bold")

        # Dimension line
        cote_y = -2.2
        ax1.annotate("", xy=(L, cote_y), xytext=(0, cote_y),
                     arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
        ax1.text(L/2, cote_y - 0.2, f"$L = {L}$ m",
                 ha="center", va="top", fontsize=10)

        for xc in [0, L]:
            ax1.plot([xc, xc], [beam_y - beam_h/2, cote_y],
                     lw=0.8, color="k", ls="--", alpha=0.4)

        # ═══════════════════════════════════════════════════════════════
        # SUBPLOT 2 — Shear
        # ═══════════════════════════════════════════════════════════════
        ax2.fill_between(x, T, 0, where=(T >= 0), alpha=0.25, color="tab:green")
        ax2.fill_between(x, T, 0, where=(T < 0), alpha=0.25, color="tab:green")
        ax2.plot(x, T, color="tab:green", lw=2.2)
        ax2.axhline(0, color="gray", lw=0.8, ls="--")

        ax2.plot(0, T_max, "o", color="tab:green", ms=6)
        ax2.annotate(f"$T_{{max}} = +{T_max:.1f}$ kN",
                    xy=(0, T_max), xytext=(0.4, T_max),
                    fontsize=9.5, color="tab:green", fontweight="bold", va="center")

        ax2.plot(L, T_min, "o", color="tab:green", ms=6)
        ax2.annotate(f"$T_{{min}} = {T_min:.1f}$ kN",
                    xy=(L, T_min), xytext=(L - 0.4, T_min),
                    fontsize=9.5, color="tab:green", fontweight="bold", va="center", ha="right")

        ax2.plot(L/2, 0, "s", color="tab:green", ms=5)
        ax2.annotate(f"$T = 0$  à  $x = {L/2:.2f}$ m",
                    xy=(L/2, 0), xytext=(L/2 + 0.6, T_max * 0.22),
                    fontsize=8.5, color="tab:green",
                    arrowprops=dict(arrowstyle="-", color="tab:green", lw=0.8))

        ax2.set_xlim(0, L)
        ax2.set_ylabel("$T(x)$  [kN]")
        ax2.set_xlabel("$x$  [m]")
        ax2.set_title("Diagramme de l'effort tranchant $T(x)$", fontweight="bold")

        # ═══════════════════════════════════════════════════════════════
        # SUBPLOT 3 — Moment
        # ═══════════════════════════════════════════════════════════════
        ax3.fill_between(x, M, 0, alpha=0.25, color="tab:purple")
        ax3.plot(x, M, color="tab:purple", lw=2.2)
        ax3.axhline(0, color="gray", lw=0.8, ls="--")

        ax3.plot(x_Mmax, M_max, "o", color="tab:purple", ms=7)
        ax3.annotate(f"$M_{{max}} = {M_max:.1f}$ kN·m\n$x = {x_Mmax:.2f}$ m",
                    xy=(x_Mmax, M_max), xytext=(x_Mmax + 3, M_max * .9),
                    fontsize=9.5, color="tab:purple", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="tab:purple", lw=1.0))

        for xc in [0, L]:
            ax3.plot(xc, 0, "s", color="tab:purple", ms=5)
        ax3.annotate("$M = 0$", xy=(0, 0), xytext=(0.4, -M_max * 0.10), fontsize=8, color="tab:purple")
        ax3.annotate("$M = 0$", xy=(L, 0), xytext=(L - 1.6, -M_max * 0.10), fontsize=8, color="tab:purple")

        ax3.set_xlim(0, L)
        ax3.set_ylabel("$M(x)$  [kN·m]")
        ax3.set_xlabel("$x$  [m]")
        ax3.set_title("Diagramme du moment fléchissant $M(x)$", fontweight="bold")

        plt.savefig(save_path, dpi=160, bbox_inches="tight")
        print(f"Saved → {save_path}")