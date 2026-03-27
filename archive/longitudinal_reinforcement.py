import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def compute_neutral_axis(M_Ed, f_cd, b_eff, h_f, b_w, h_tot, d):
    # https://moodle.insa-lyon.fr/pluginfile.php/438753/mod_resource/content/1/5%20-%20Flexion%20sans%20charge%20axiale%20ELU.pdf
    # slide 28
    M_ct = b_eff * h_f * f_cd * (d - h_f / 2)  # N·m
    if M_Ed < M_ct:
        # N.A. in the flange (table)
        print("N.A. in the flange (table)")
        # Design a rectangular section (beff x h) balanced MEd
        # Slide 25 : https://moodle.insa-lyon.fr/pluginfile.php/438753/mod_resource/content/1/5%20-%20Flexion%20sans%20charge%20axiale%20ELU.pdf
        mu = M_Ed / (b_eff * d**2 * f_cd)  # Reduced moment μ dimensionless
        alpha = (1 - np.sqrt(1 - 2 * mu)) / 0.8
        x_na = alpha * d
        if mu > 0.372 or alpha > 0.617:
            print("Resize the concrete section (b x h) or increase fck or add A’s")
    else:
        # N.A. in the web (âme)
        print("N.A. in the web (âme)")
        raise (ValueError("Not implemented, not useful in our case"))
    return dict(M_ct=M_ct, mu=mu, alpha=alpha, x_na=x_na)


def calc_sigma_s(f_yd, epsilon_s, Es=200e9):
    """
    Steel stress-strain law for HA B500B (inclined top branch, EC2).

    Parameters
    ----------
    epsilon_s : float — steel strain [dimensionless, e.g. 0.00217]
    Es        : float — Young's modulus [Pa], default 200 GPa

    Returns
    -------
    sigma : float — steel stress [Pa]

    The law has two branches (slide stress-strain diagram):
      - Elastic  : ε ≤ ε_se  →  σ = Es · ε          (linear)
      - Hardening: ε > ε_se  →  σ = 435e6 + 0.817e6·(ε·1000 - 2.17)  (inclined plateau)
    ε_se = f_yd / Es = 435e6 / 200e9 ≈ 2.17 ‰ is the yield strain.
    The hardening slope (0.817e6 Pa/‰) models the inclined top branch of HA B500B.
    """
    epsilon_se = f_yd / Es  # yield strain [dimensionless]

    if epsilon_s <= epsilon_se:
        return Es * epsilon_s
    else:
        # epsilon_s * 1000 converts to ‰ to match the 0.817e6 [Pa/‰] slope
        return f_yd + 0.817e6 * (epsilon_s * 1000 - 2.17)


def plot_neutral_axis(
    M_Ed,
    f_cd,
    f_ctm,
    f_yk,
    f_yd,
    h_tot,
    b_eff,
    h_f,
    b_w,
    d,
    save_path="./neutral_axis.svg",
):
    """
    # Tracer le profil en fonction de z sigma_x = Mz/EIz y
    # https://moodle.insa-lyon.fr/pluginfile.php/438753/mod_resource/content/1/5%20-%20Flexion%20sans%20charge%20axiale%20ELU.pdf
    # Page 16
    #"""

    h_w = h_tot - h_f
    res = compute_neutral_axis(M_Ed, f_cd, b_eff, h_f, b_w, h_tot, d)
    print(res)
    x_na = res["x_na"]
    alpha = res["alpha"]
    y_na = h_tot - x_na  # y of N.A. from bottom
    y_as = h_tot - d  # y of steel from bottom

    epsilon_c = 3.5e-3
    epsilon_s = epsilon_c * (1 - alpha) / alpha  # steel strain (tension)
    sigma_s = calc_sigma_s(f_yd, epsilon_s)
    print(epsilon_c, epsilon_s, sigma_s)

    F_cf = (b_eff - b_w) * h_f * f_cd  # N  — flange overhangs
    F_cw = 0.8 * b_w * alpha * d * f_cd  # N  — web

    if M_Ed <= res["M_ct"]:
        # N.A. in flange → single rectangular block b_eff (slide 20)
        F_c = 0.8 * alpha * b_eff * d * f_cd
        A_s = F_c / sigma_s
    else:
        # N.A. in web → T decomposition (slide 32)
        A_s = (F_cf + F_cw) / sigma_s

    A_s_min = max(0.26 * f_ctm / f_yk * b_w * d, 0.0013 * b_w * d)

    # Check As,min
    # Select bars diameter
    # Check dreal ≥ d and d'real ≤ d'
    print(f"A_s = {A_s*1e4:.1f} cm²")
    print(f"A_s_min = {A_s_min*1e4:.1f} cm²")
    if A_s < A_s_min:
        print("A_s < A_s_min")
    else:
        print("check : A_s >= A_s_min")

    # ε(y) linear: zero at y_na, +epsilon_c at top (y=h_tot), -epsilon_s at bottom (y=y_as)
    # ε(y) = epsilon_c * (y - y_na) / x_na
    y_profile = np.array([0, y_as, y_na, h_tot])
    eps_profile = epsilon_c * (y_profile - y_na) / x_na

    fig, (ax_sec, ax_eps) = plt.subplots(
        1,
        2,
        sharey=True,
        figsize=(8, 6),
        gridspec_kw={"width_ratios": [3, 1.5]},
    )
    fig.suptitle(f"Position de l'axe neutre : y = {x_na*100:.1f} cm", fontsize=11)

    # ── Section ───────────────────────────────────────────────────────────
    # ax_sec.set_aspect("equal", adjustable="datalim")
    # ax_sec.set_aspect("equal")
    ax_sec.add_patch(patches.Rectangle((-b_eff / 2, h_w), b_eff, h_f))
    ax_sec.add_patch(patches.Rectangle((-b_w / 2, 0), b_w, h_w))
    ax_sec.axhline(
        y_na, color="k", lw=1.5, ls="--", label=f"A.N.  y = {x_na*100:.1f} cm"
    )
    for xb in np.linspace(-b_w / 2 + 0.06, b_w / 2 - 0.06, 6):
        ax_sec.add_patch(
            plt.Circle((xb, y_as), 0.018, fc="k", zorder=5)  # pyright: ignore[reportPrivateImportUsage]
        )  # pyright: ignore
    ax_sec.plot([], [], "ko", ms=6, label="Armatures As")
    ax_sec.set_xlabel("z [m]")
    ax_sec.set_ylabel("y [m]")
    ax_sec.legend(fontsize=8)

    # Synchroniser l'axe y du diagramme de déformation avec la section
    ax_eps.set_ylim(ax_sec.get_ylim())
    ax_eps.set_yticks(ax_sec.get_yticks())
    ax_eps.set_yticklabels([f"{tick:.2f}" for tick in ax_sec.get_yticks()])

    # ── Strain ────────────────────────────────────────────────────────────
    ax_eps.plot(eps_profile * 1000, y_profile, "k-", lw=1.5)
    ax_eps.axvline(0, color="k", lw=0.8, ls=":")
    ax_eps.axhline(y_na, color="k", lw=1.5, ls="--")
    ax_eps.set_xlabel("ε [‰]")

    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Saved → {save_path}")
