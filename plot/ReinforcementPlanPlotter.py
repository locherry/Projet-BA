from __future__ import annotations
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Arc
import math
import numpy as np


def _ensure_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


class ReinforcementPlanPlotter:

    def __init__(
        self,
        curtailment,
        design,
        beam,
        shear,
        save_path: str = "./plot/reinforcement_plan.svg",
    ):
        self.curtailment = curtailment
        self.design = design
        self.beam = beam
        self.shear = shear          # ShearDesign instance — optional
        self.save_path = save_path

    @property
    def L(self) -> float:
        return self.beam.L

    @property
    def section(self):
        return self.design.section

    # ─────────────────────────────────────────────────────────────

    def plot(self):
        L = self.L
        section = self.section

        h_tot = section.h_tot   # total height [m]
        h_f   = section.h_f     # flange thickness [m]
        h_w   = h_tot - h_f     # web height [m]
        b_eff = section.b_eff   # effective flange width [m]
        b_w   = section.b_w     # web width [m]

        # ── Axes limits ───────────────────────────────────────────
        # y: 0 = bottom fibre, h_tot = top fibre (real metres)
        # x: 0 → L  (real metres)
        margin_top    = h_tot * 0.35   # room for global length arrow + title
        margin_bottom = h_tot * 0.55   # room for dimension lines below beam
        
        cotation_offset = 0.4 * h_tot

        fig, ax = plt.subplots(figsize=(16, 5))

        ax.set_xlim(0, L)
        ax.set_ylim(-margin_bottom, h_tot + margin_top)
        ax.set_aspect("auto")           # keep x/y independently scaled
        ax.set_xlabel("x [m]")
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.6)

        # ── Draw T-section beam outline (longitudinal view) ────────
        # The beam is shown as a filled rectangle: full width = b_eff (flange)
        # at top, narrower web = b_w at bottom, centred on x-axis.
        # Since this is a plan view we draw the elevation (side profile):
        # just the web rectangle + flange rectangle stacked.

        # Flange: spans full beam length, height h_f at top
        flange_y_bottom = h_w  # from bottom fibre

        # We draw the beam silhouette as two overlaid patches:
        # — web (b_w wide — represented symbolically as the full length elevation)
        # — flange (b_eff wide — shown via slightly darker top band)

        # Web zone (full length, bottom → h_w)
        ax.add_patch(
            patches.Rectangle(
                (0, 0), L, h_w,
                fc="0.93", ec="0.35", lw=1.0, zorder=0,
            )
        )
        # Flange zone (full length, h_w → h_tot)
        ax.add_patch(
            patches.Rectangle(
                (0, h_w), L, h_f,
                fc="0.85", ec="0.35", lw=1.0, zorder=0,
            )
        )

        # Flange/web interface dashed line
        ax.axhline(h_w, color="0.5", linestyle="--", linewidth=0.6, zorder=1)
        
        
        # ── Draw stirrups ─────────────────────────────────────────
        self._draw_stirrups(ax, h_w, h_tot, cotation_offset)

        # ── Draw reinforcement bars ────────────────────────────────
        # Each layer is plotted at its actual y position (from bottom fibre),
        # taken from the bar positions in the layout.

        # Build a mapping: layer index → y coordinate of bar centroid
        bar_positions = self.design.layout.get_bar_positions()

        # Group bar_positions by row (they are ordered bottom → top within rows)
        # We need one representative y per curtailment layer.
        # Curtailment layers are defined per row in the layout.
        # Simplest: use the y of the first bar in each layout row.

        # Collect unique y values (one per distinct diameter / row)
        # bar_positions: list of (x_cross, y_from_bottom, radius)
        # We want unique y levels grouped by radius (diameter)
        seen = {}
        y_by_diameter = {}  # phi -> y (from bottom)
        for _, y_bar, r in bar_positions:
            phi = round(r * 2, 6)
            if phi not in y_by_diameter:
                y_by_diameter[phi] = y_bar

        # Map curtailment layers to y positions
        layer_y = []
        for layer in self.curtailment.layers:
            phi = round(layer.phi, 6)
            y_bar = y_by_diameter.get(phi, h_w / 2)  # fallback to mid-web
            layer_y.append(y_bar)

        for (layer, cutoff), y in zip(
            zip(self.curtailment.layers, self.curtailment.cutoffs),
            layer_y,
        ):
            phi = layer.phi
            r_hook = 4 * phi
            v_hook = 5 * phi
            if phi == 40e-3:
                # v_hook = 982e-3
                r_hook = 200e-3
            elif phi == 32e-3:
                # v_hook = 646e-3
                r_hook = 145e-3
                


            # Bar extents (same logic as original)
            xL = cutoff.x_left  + r_hook
            xR = cutoff.x_right - r_hook
            if phi == 40e-3:
                xL = 0 + r_hook + self.design.layout.c_nom
                xR = L - r_hook - self.design.layout.c_nom

            # Color mapping
            if "40" in layer.label:
                color = "tab:blue"
            elif "32" in layer.label:
                color = "tab:orange"
            else:
                color = "tab:gray"

            lw = 2.5

            # Main bar line
            ax.plot([xL, xR], [y, y], linewidth=lw, color=color, zorder=3)

            # Curtailment cut-off vertical dashed lines
            ax.axvline(cutoff.x_left,  linestyle="--", linewidth=0.8,
                       color="0.4", zorder=2)
            ax.axvline(cutoff.x_right, linestyle="--", linewidth=0.8,
                       color="0.4", zorder=2)

            # Hooks
            self._draw_hook(ax, xL, y, "left",  r_hook, v_hook, lw, color)
            self._draw_hook(ax, xR, y, "right", r_hook, v_hook, lw, color)

            # Bar label (right of bar)
            ax.text(
                xR + r_hook + 0.005 * L, y,
                layer.label,
                va="center", fontsize=7.5, color=color, fontweight="bold",
                zorder=4,
            )

            # ── Bar length dimension ──────────────────────────────
            y_dim = y - 0.04 * h_tot
            bar_left  = xL - r_hook
            bar_right = xR + r_hook
            ax.annotate(
                "",
                xy=(bar_right, y_dim - 0.8*cotation_offset), xytext=(bar_left, y_dim - 0.8*cotation_offset),
                arrowprops=dict(arrowstyle="<->", linewidth=0.7),
                zorder=4,
            )
            ax.text(
                (bar_left + bar_right) / 2, y_dim - 0.025 * h_tot - 0.9*cotation_offset,
                f"{bar_right - bar_left:.2f} m",
                ha="center", fontsize=6.5, zorder=4,
            )

            # ── Anchorage lengths ─────────────────────────────────
            y_lbd = y - 0.09 * h_tot

            # left anchorage
            ax.annotate(
                "",
                xy=(bar_left, y_lbd - cotation_offset),
                xytext=(bar_left + cutoff.lbd, y_lbd - cotation_offset),
                arrowprops=dict(arrowstyle="<->", linewidth=0.7, color="tab:red"),
                zorder=4,
            )
            ax.text(
                bar_left + cutoff.lbd / 2, y_lbd - 0.025 * h_tot - 1.10*cotation_offset,
                f"l_bd={cutoff.lbd*100:.0f} cm",
                ha="center", fontsize=6, color="tab:red", zorder=4,
            )

            # right anchorage
            ax.annotate(
                "",
                xy=(bar_right, y_lbd - cotation_offset),
                xytext=(bar_right - cutoff.lbd, y_lbd - cotation_offset),
                arrowprops=dict(arrowstyle="<->", linewidth=0.7, color="tab:red"),
                zorder=4,
            )
            ax.text(
                bar_right - cutoff.lbd / 2, y_lbd - 0.025 * h_tot - 1.10*cotation_offset,
                f"l_bd={cutoff.lbd*100:.0f} cm",
                ha="center", fontsize=6, color="tab:red", zorder=4,
            )

        # ── Global beam length arrow ───────────────────────────────
        y_top = h_tot + 0.1 * h_tot
        ax.annotate(
            "",
            xy=(L, y_top), xytext=(0, y_top),
            arrowprops=dict(arrowstyle="<->", linewidth=0.9),
        )
        ax.text(
            L / 2, y_top + 0.04 * h_tot,
            f"L = {L:.2f} m",
            ha="center", fontsize=9, fontweight="bold",
        )

        # ── y-axis labels (heights) ────────────────────────────────
        ax.axhline(0,     color="0.3", linewidth=0.8)   # bottom fibre
        ax.axhline(h_tot, color="0.3", linewidth=0.8)   # top fibre

        ax.text(-0.01 * L, 0,     "±0",        va="center", ha="right", fontsize=7, color="0.4")
        ax.text(-0.01 * L, h_tot, f"{h_tot*100:.0f} cm", va="center", ha="right", fontsize=7, color="0.4")
        ax.text(-0.01 * L, h_w,   f"{h_w*100:.0f} cm",  va="center", ha="right", fontsize=7, color="0.4")

        ax.set_title("Plan de ferraillage longitudinal", fontsize=12, fontweight="bold", pad=12)

        plt.tight_layout()
        _ensure_dir(self.save_path)
        plt.savefig(self.save_path, dpi=150, bbox_inches="tight")
        print(f"Saved → {self.save_path}")
        plt.close(fig)

    # ─────────────────────────────────────────────────────────────

    def _draw_stirrups(self, ax, h_w: float, h_tot: float, cotation_offset: float):
        shear  = self.shear
        phi_t  = shear.phi_t
        L      = self.L
        c_nom  = self.design.layout.c_nom
        y_bot  = c_nom
        y_top  = h_tot - c_nom
        color  = "tab:green"
        lw     = 1.4

        all_zones = shear.stirrup_zones()
        zones_left = [z for z in all_zones if z["x_end"] <= L / 2 + 1e-9]

        def draw_stirrup(x_draw):
            ax.plot([x_draw, x_draw], [y_bot, y_top],
                    color=color, linewidth=lw, zorder=2)
            ax.plot([x_draw - phi_t/2, x_draw + phi_t/2], [y_bot, y_bot],
                    color=color, linewidth=lw, zorder=2)
            ax.plot([x_draw - phi_t/2, x_draw + phi_t/2], [y_top, y_top],
                    color=color, linewidth=lw, zorder=2)
        
        ax.plot([0, L], [y_top, y_top],
                color=color, linewidth=lw, zorder=2)

        # Collect all left-half stirrup x positions, building from left edge inward
        xs_left = []  # list of (x, s) to track spacing per position
        x_cursor = 0.0

        for zone in zones_left:
            s       = zone["s"]
            x_start = zone["x_start"]
            x_end   = zone["x_end"]

            # First zone starts at 0, subsequent zones start right after last placed
            if not xs_left:
                x_cursor = x_start  # = 0 for first zone
            # else: x_cursor is already just past the previous zone's last stirrup

            while x_cursor <= x_end + 1e-9:
                xs_left.append((x_cursor, s))
                x_cursor += s

            # Zone boundary + annotation
            if x_start > 0:
                ax.axvline(x_start, color=color, linewidth=0.5,
                        linestyle=":", alpha=0.5, zorder=1)
                ax.axvline(L - x_start, color=color, linewidth=0.5,
                        linestyle=":", alpha=0.5, zorder=1)

            # Alternate annotation level to avoid overlaps
            level = zones_left.index(zone) % 2
            y_ann = y_bot - 0.06 * h_tot - (0.3 + 0.25 * level) * cotation_offset
            x_mid_left  = (x_start + x_end) / 2
            x_mid_right = L - x_mid_left
            for x_mid in [x_mid_left, x_mid_right]:
                ax.text(x_mid, y_ann, f"st={s*100:.0f} cm",
                        ha="center", fontsize=6, color=color, zorder=4)
        # Draw all stirrups: left positions + their exact mirrors
        # The center gap (between last left stirrup and its mirror) is filled
        # naturally — no overlap because mirror of xs > L/2 when xs < L/2
        drawn = set()
        for xs, s in xs_left:
            for x_draw in [xs, L - xs]:
                key = round(x_draw, 6)
                if key not in drawn:
                    draw_stirrup(x_draw)
                    drawn.add(key)
                    
        draw_stirrup(L / 2)
                    
        # Global label
        ax.text(
            L / 2, y_bot - 0.06 * h_tot - 0.35 * cotation_offset,
            f"ø{phi_t*1e3:.0f} mm",
            ha="center", fontsize=6.5, color=color, fontweight="bold", zorder=4,
        )
        
    def _draw_hook(self, ax, x, y, side, r, v, lw, color):
        """Draw a 45° hook tangent to the bar end."""

        if side == "left":
            xc = x + r

            arc = Arc(
                (xc - r, y + r),
                2 * r,
                2 * r,
                theta1=180 - 45,
                theta2=270,
                linewidth=lw,
                color=color,
            )
            ax.add_patch(arc)

            # ── extension line (45° upward-right)
            angle = math.radians(45)

            cx = x    # circle center x
            cy = y + r         # circle center y

            # exact arc endpoint at theta2 = 270°
            x_end = cx - r* math.cos(math.radians(45))
            y_end = cy + r* math.cos(math.radians(45))

            x2 = x_end + v * math.cos(angle)
            y2 = y_end + v * math.sin(angle)

            ax.plot([x_end, x2], [y_end, y2], linewidth=lw, color=color)


        elif side == "right":
            xc = x - r

            arc = Arc(
                (xc + r, y + r),
                2 * r,
                2 * r,
                theta1=270,
                theta2=360 + 45,
                linewidth=lw,
                color=color,
            )
            ax.add_patch(arc)

            # ── extension line (45° upward-right)
            angle = math.radians(45+90)

            cx = x    # circle center x
            cy = y + r         # circle center y

            # exact arc endpoint at theta2 = 270°
            x_end = cx + r* math.cos(math.radians(45))
            y_end = cy + r* math.cos(math.radians(45))

            x2 = x_end + v * math.cos(angle)
            y2 = y_end + v * math.sin(angle)

            ax.plot([x_end, x2], [y_end, y2], linewidth=lw, color=color)