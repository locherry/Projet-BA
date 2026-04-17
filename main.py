from classes.TSection import TSection
from classes.Materials import Concrete, Steel
from classes.Beam import Beam, DistributedLoad
from classes.FlexuralDesign import FlexuralDesign
from classes.LoadCalculator import LoadCalculator
from classes.ReinforcementLayout import ReinforcementLayout
from classes.Rules import DurabilityRules

from plot.BeamPlotter import BeamPlotter
from plot.SectionPlotter import SectionPlotter
import numpy as np

# ---------------------------------------------------------------------------- #
#                                 Define Params                                #
# ---------------------------------------------------------------------------- #

# Loads
load_calc = LoadCalculator()
q_elu = load_calc.q_ELU          # ULS distributed load [N/m]
q_eqp = load_calc.q_EQP()        # quasi-permanent SLS load [N/m]

# Geometry
L = 15.25  # [m]

beam_uls = Beam(L, DistributedLoad(q_elu))
beam_eqp = Beam(L, DistributedLoad(q_eqp))

section = TSection(
    h_tot=1.25,
    b_eff=2.0,
    h_f=0.15,
    b_w=0.5,
)

concrete = Concrete(f_ck=30e6, gamma_c=1.5)
steel    = Steel(f_yk=500e6,  gamma_s=1.15)

# First estimate of effective depth — will be updated to d_reel once
# the reinforcement layout is known.
d_estimate = 0.9 * section.h_tot

design = FlexuralDesign(section, concrete, steel, d=d_estimate)

M_Ed  = beam_uls.max_moment()
M_Eqp = beam_eqp.max_moment()

# Durability parameters
rules = DurabilityRules()
work_life_years = 50
exposure_class = "XC1"
c_min_dur = rules.c_min_dur(exposure_class, work_life_years)
w_max = rules.w_max(exposure_class)

# ---------------------------------------------------------------------------- #
#                       Step 1 — Preliminary A_s (with d estimate)            #
# ---------------------------------------------------------------------------- #

A_s_required = design.reinforcement_area(M_Ed)
print(f"M_Ed        = {M_Ed/1e6:.2f} MN·m")
print(f"M_Eqp       = {M_Eqp/1e6:.2f} MN·m")
print(f"A_s required (d estimate) = {A_s_required*1e4:.2f} cm²")

# ---------------------------------------------------------------------------- #
#                       Step 2 — Define layout and update d                    #
# ---------------------------------------------------------------------------- #

# Stirrup diameter is a design choice — ø10 mm here.
layout = ReinforcementLayout(section, phi_t=10e-3)

layout.set_exposure(c_min_dur)

# Rows declared bottom → top
# Bottom row: 4 HA 32
layout.add_row(n_groups=4, bars_per_group=1, diameter=32e-3, grouped=True)
# Top row: 4 HA 25
layout.add_row(n_groups=4, bars_per_group=1, diameter=25e-3, grouped=False)

# Attach layout — this also updates design.d to d_reel
design.set_layout(layout)

d_reel = design.d
print(f"\nd_reel      = {d_reel*1e2:.2f} cm  (replaces d estimate)")

# ---------------------------------------------------------------------------- #
#                       Step 3 — Recompute A_s with d_reel                    #
# ---------------------------------------------------------------------------- #

A_s = design.reinforcement_area(M_Ed)
print(f"A_s required (d_reel)     = {A_s*1e4:.2f} cm²")

# Provided area (4 × ø32 bottom + 4 × ø25 top)
A_s_provided = 4 * np.pi * (32e-3 / 2) ** 2 + 4 * np.pi * (25e-3 / 2) ** 2
print(f"A_s provided              = {A_s_provided*1e4:.2f} cm²")

# ---------------------------------------------------------------------------- #
#                       Step 4 — SLS Stress Limitation                        #
# ---------------------------------------------------------------------------- #

stress = design.check_stress_limitation(M_Ed, M_Eqp, A_s=A_s_provided)

print("\n--- SLS Stress Limitation (EC2 §7.2) ---")
print(f"  φ(∞,t0)   = {stress['phi_inf_t0']:.2f}")
print(f"  φ_ef      = {stress['phi_ef']:.2f}")
print(f"  E_c,eff   = {stress['E_c_eff']/1e9:.2f} GPa")
print(f"  n (αe)    = {stress['n']:.2f}")
print(f"  x_SLS     = {stress['x_SLS']*1e2:.2f} cm")
print(f"  σ_c       = {stress['sigma_c']/1e6:.2f} MPa  (limit {stress['limit_compression']/1e6:.1f} MPa)  "
      f"{'✓' if stress['ok_compression'] else '✗'}")
print(f"  σ_s       = {stress['sigma_s']/1e6:.2f} MPa  (limit {stress['limit_steel']/1e6:.1f} MPa)  "
      f"{'✓' if stress['ok_steel'] else '✗'}")

# ---------------------------------------------------------------------------- #
#                       Step 5 — Crack Control                                 #
# ---------------------------------------------------------------------------- #
crack = design.crack_control(M_Ed, M_Eqp, A_s=A_s_provided, w_max=w_max)

print("\n--- Crack Control (EC2 §7.3) ---")
print(f"  As,min    = {crack['As_min']*1e4:.2f} cm²  {'✓' if crack['ok_As_min'] else '✗'}")
print(f"  x (SLS)   = {crack['x_SLS']*1e2:.2f} cm")
print(f"  σ_s       = {crack['sigma_s']/1e6:.1f} MPa")
print(f"  h_c,ef    = {crack['h_c_ef']*1e2:.2f} cm")
print(f"  ρ_p,eff   = {crack['rho_p_eff']:.4f}")
print(f"  φ_eq      = {crack['phi_eq']*1e3:.1f} mm")
print(f"  s_r,max   = {crack['s_r_max']*1e3:.1f} mm")
print(f"  w_k       = {crack['w_k']*1e3:.3f} mm  (limit {crack['w_max']*1e3:.1f} mm)  "
      f"{'✓' if crack['ok_wk'] else '✗'}")

# ---------------------------------------------------------------------------- #
#                                   Plots                                      #
# ---------------------------------------------------------------------------- #

beam_plotter = BeamPlotter(beam_uls)
beam_plotter.plot()

section_plotter = SectionPlotter(design)
section_plotter.plot(M_Ed)