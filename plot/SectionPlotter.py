import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class SectionPlotter:
    def __init__(self, design):
        self.design = design  # FlexuralDesign

    def plot(self, M_Ed, save_path="./plot/section.svg"):
        section = self.design.section

        # ── Get results from design (NO recomputation)
        res = self.design.neutral_axis(M_Ed)

        x_na = res["x_na"]

        h_tot = section.h_tot
        y_na = h_tot - x_na
        y_as = h_tot - self.design.d

        # ── Strain computation (EC2)
        epsilon_c = 3.5e-3

        # ── Reinforcement (reuse design method)
        A_s = self.design.reinforcement_area(M_Ed)

        # ── Strain profile
        y_profile = np.array([0, y_as, y_na, h_tot])
        eps_profile = epsilon_c * (y_profile - y_na) / x_na

        # ── Plot
        fig, (ax_sec, ax_eps) = plt.subplots(
            1,
            2,
            sharey=True,
            figsize=(8, 6),
            gridspec_kw={"width_ratios": [3, 1.5]},
        )

        fig.suptitle(
            f"Axe neutre : x = {x_na*100:.1f} cm | A_s = {A_s*1e4:.1f} cm²", fontsize=11
        )

        # ═══════════════════════════════════════════════════════════════
        # SECTION
        # ═══════════════════════════════════════════════════════════════
        h_w = section.h_w

        ax_sec.add_patch(
            patches.Rectangle((-section.b_eff / 2, h_w), section.b_eff, section.h_f)
        )
        ax_sec.add_patch(patches.Rectangle((-section.b_w / 2, 0), section.b_w, h_w))

        # Neutral axis
        ax_sec.axhline(
            y_na, color="k", lw=1.5, ls="--", label=f"A.N. x = {y_na:.2f} cm"
        )

        layout = self.design.layout  # assuming you store it here

        bars = layout.get_bar_positions()

        for x, y, r in bars:
            ax_sec.add_patch(
                plt.Circle( # pyright: ignore[reportPrivateImportUsage]
                    (x, y), r, fc="k", zorder=5
                )
            )

        ax_sec.plot([], [], "ko", ms=6, label="Armatures As")
        # c_nom = layout.c_nom
        # phi_t = layout.phi_t
        # rect = patches.Rectangle(
        #     (-section.b_w / 2 + c_nom, c_nom),
        #     section.b_w - 2 * c_nom,
        #     section.h_tot - 2 * c_nom,
        #     fill=False,
        #     linestyle="--",
        #     linewidth=phi_t * 100,  # convert meters → points for visibility
        # )
        # ax_sec.add_patch(rect)

        ax_sec.set_xlabel("z [m]")
        ax_sec.set_ylabel("y [m]")
        ax_sec.legend(fontsize=8)

        # ═══════════════════════════════════════════════════════════════
        # STRAIN DIAGRAM
        # ═══════════════════════════════════════════════════════════════
        ax_eps.plot(eps_profile * 1000, y_profile, "k-", lw=1.5)
        ax_eps.axvline(0, color="k", lw=0.8, ls=":")
        ax_eps.axhline(y_na, color="k", lw=1.5, ls="--")

        ax_eps.set_xlabel("ε [‰]")

        # Sync y-axis
        ax_eps.set_ylim(ax_sec.get_ylim())
        ax_eps.set_yticks(ax_sec.get_yticks())
        ax_eps.set_yticklabels([f"{tick:.2f}" for tick in ax_sec.get_yticks()])

        plt.tight_layout()
        plt.savefig(save_path)
        print(f"Saved → {save_path}")
