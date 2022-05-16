import matplotlib.pyplot as plt
from libdenavit.section import FiberSection, FiberCirclePatch, FiberSingle, AciStrainCompatibility
import numpy as np
from math import cos, sin, pi

# Units: kips, in, ksi
D  = 10     # Diameter
fc = 4      # Concrete compressive strength
fy = 60     # Steel yield strength
Es = 29000  # Steel modulus of elasticity 

# Build FiberSection Object
c = FiberSection(100, 100)
b1 = FiberCirclePatch(0, 0, 0, D / 2, 1)
a1 = FiberSingle(.655, 0, -.3 * D, 2, 1)  # A, x, y, m, m_neg
a2 = FiberSingle(.655, 0, .3 * D, 2, 1)
a3 = FiberSingle(.655, sin(pi / 3) * .3 * D, -cos(pi / 3) * .3 * D, 2, 1)
a4 = FiberSingle(.655, sin(pi / 3) * .3 * D, cos(pi / 3) * .3 * D, 2, 1)
a5 = FiberSingle(.655, -sin(pi / 3) * .3 * D, -cos(pi / 3) * .3 * D, 2, 1)
a6 = FiberSingle(.655, -sin(pi / 3) * .3 * D, cos(pi / 3) * .3 * D, 2, 1)
c.add_fibers(b1, a1, a2, a3, a4, a5, a6)

# Output Information From FiberSection
c.print_section_properties()
c.plot_fibers(3)

# Build ACI Strain Compatibility Object
test = AciStrainCompatibility(c)
test.add_concrete_boundary(0.5 * D, 0, 0)
test.add_concrete_boundary(-0.5 * D, 0, 0)
test.add_concrete_boundary(0, -0.5 * D, 0)
test.add_concrete_boundary(0, 0.5 * D, 0)
test.add_steel_boundary(3, 3, 0)
test.add_steel_boundary(3, -3, 0)
test.add_steel_boundary(-3, 3, 0)
test.add_steel_boundary(-3, -3, 0)
test.add_material(1, "concrete", fc, "US")
test.add_material(2, "steel", fy, Es)
test.build_data()

# Select Neutral Axis Locations
a = np.arange(-70, -30, 0.5)
b = np.arange(-30, 20, 0.1)
c = np.arange(20, 80, 0.1)
d = np.concatenate((a, b, c))

# Perform Interaction Calculations
plot_P = []
plot_Mx = []
for i in d:
    P, Mx, My, et = test.compute_point(0, i, 0)
    plot_P.append(-P)
    plot_Mx.append(Mx)

# Plot Interaction Diagram
plt.plot(plot_Mx, plot_P, 'bo-')
plt.show()