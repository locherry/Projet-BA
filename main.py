from classes.TSection import TSection
from classes.Materials import Concrete, Steel
from classes.Beam import Beam, DistributedLoad
from classes.FlexuralDesign import FlexuralDesign
from classes.LoadCalculator import LoadCalculator
from classes.ReinforcementLayout import ReinforcementLayout

from plot.BeamPlotter import BeamPlotter
from plot.SectionPlotter import SectionPlotter
import numpy as np

# Load
load_calc = LoadCalculator()
q = load_calc.q_ELU # N/m

# Geometry
L = 15.25

beam = Beam(L, DistributedLoad(q))

section = TSection(
    h_tot=1.25,
    b_eff=2.0,
    h_f=0.15,
    b_w=0.5
)

concrete = Concrete(f_ck=30e6, gamma_c=1.5)
steel = Steel(f_yk=500e6, gamma_s=1.15)

d = 0.9 * section.h_tot

design = FlexuralDesign(section, concrete, steel, d)

M_Ed = beam.max_moment()

A_s = design.reinforcement_area(M_Ed)

print(f"M_Ed = {M_Ed/1e6:.2f} MN·m")
print(f"A_s = {A_s*1e4:.2f} cm²")

layout = ReinforcementLayout(section)
design.layout = layout

# Bottom row: 3 groups of 2 bars
layout.add_row(n_groups=3, bars_per_group=2, diameter=32e-3, grouped=True)

# Top row: 1 bar
layout.add_row(n_groups=1, bars_per_group=1, diameter=32e-3, grouped=False)

d_reel = layout.compute_d_reel()
print(f"d_reel = {d_reel*1e2:.2f} cm")

beam_plotter = BeamPlotter(beam)
beam_plotter.plot()

section_plotter = SectionPlotter(design)
section_plotter.plot(M_Ed)