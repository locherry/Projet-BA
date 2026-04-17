import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class ReinforcementPlotter:
    def __init__(self, design):
        self.design = design  # FlexuralDesign


    def plot(self, M_Ed, save_path="./plot/reinforcement.svg"):
        section = self.design.section
        layout  = self.design.layout   # ReinforcementLayout


        # ── Geometry
        h_tot = section.h_tot
        h_w   = section.h_w
        b_eff = section.b_eff
        b_w   = section.b_w
        c_nom = layout.c_nom
        phi_t = layout.phi_t


        # ── Figure
        fig, ax_sec = plt.subplots(figsize=(7, 6))
        ax_sec.set_xlabel("z [m]")
        ax_sec.set_ylabel("y [m]")
        ax_sec.set_title("Armatures longitudinales (vue de section)", fontsize=11)


        # ── Concrete
        # flange
        ax_sec.add_patch(
            patches.Rectangle(
                (-b_eff / 2, h_w),
                b_eff,
                section.h_f,
                fc="0.9",
                ec="0.3",
                lw=0.8,
            )
        )
        # web
        ax_sec.add_patch(
            patches.Rectangle(
                (-b_w / 2, 0),
                b_w,
                h_w,
                fc="0.9",
                ec="0.3",
                lw=0.8,
            )
        )


        # ── Main bars
        bars = layout.get_bar_positions()

        for x, y, r in bars:
            ax_sec.add_patch(
                patches.Circle((x, y), r, fc="k", ec="none", zorder=5)
            )

        ax_sec.plot([], [], "ko", ms=6, label="Armatures principales As")


        # ── 1) External stirrup (perimeter)
        stir_zmin = -b_w / 2 + c_nom
        stir_zmax =  b_w / 2 - c_nom
        stir_ymin =  c_nom + phi_t
        stir_ymax =  h_tot - c_nom

        ax_sec.add_patch(
            patches.Rectangle(
                (stir_zmin, stir_ymin),
                stir_zmax - stir_zmin,
                stir_ymax - stir_ymin,
                fc="none",
                ec="0.4",
                lw=1.2,
            )
        )


        # ── 2) Top‑bar level and top stirrup bar

        if not bars:
            raise ValueError("No bars in ReinforcementLayout")

        # top‑bar level
        y_top = max(y for x, y, r in bars)+ phi_t
        tol   = 1e-3
        top_bars = [b for b in bars if abs(b[1] - y_top) < tol]

        # horizontal stirrup bar at y_top, fixed to outer stirrup legs
        ax_sec.plot(
            [stir_zmin, stir_zmax],
            [y_top, y_top],
            color="0.25",
            lw=1.2,
            zorder=4,
        )


        # ── Legend
        ax_sec.plot([], [], "-", lw=1.2, c="0.4",  label="Cadre extérieur (périmétrique)")
        ax_sec.plot([], [], "-", lw=1.2, c="0.25", label="Barre horizontale de maintien (top)")

        ax_sec.legend()
        plt.tight_layout()
        plt.savefig(save_path)
        print(f"Saved → {save_path}")