import matplotlib.pyplot as plt
from libdenavit.section import RC,Obround,ReinfIntersectingLoops
import numpy as np

# Units: kips, in, ksi
D  = 48     # Cross-sectional width
a  = 24     # Cross-sectional height
fc = 4      # Concrete compressive strength
fy = 60     # Steel yield strength
Es = 29000
nb = 36
Ab = 1.0
axis = 'x'
units = 'us'

# Define RC section object
conc_cross_section = Obround(D, a)
reinforcement = ReinfIntersectingLoops(D-6, a, 36, 1.0)
section = RC(conc_cross_section, reinforcement, fc, fy, 'US')

# Plot Section
#section.plot_section(show=True)

aci_obj = section.aci_strain_compatibility_object()
aci_obj.fiber_section.plot_fibers()
aci_obj.fiber_section.print_section_properties()

print(f'A  = {conc_cross_section.A}')
print(f'Ix = {conc_cross_section.Ix}')
print(f'Iy = {conc_cross_section.Iy}')