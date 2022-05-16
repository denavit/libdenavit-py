import matplotlib.pyplot as plt
from libdenavit.section import FiberSection, FiberQuadPatch, FiberSingle, AciStrainCompatibility
import numpy as np

# Units: kips, in, ksi
B  = 20     # Cross-sectional width
H  = 40     # Cross-sectional height
fc = 4      # Concrete compressive strength
fy = 60     # Steel yield strength
Es = 29000  # Steel modulus of elasticity 

# Build Fiber Section Object
c = FiberSection(50, 50)
b1 = FiberQuadPatch(-10, -20, -10, 20, 10, 20, 10, -20, 1)  # xI, yI, ..., m
a1 = FiberSingle(12, -4, -14, 2, 1)  # A, x, y, m, m_neg
a2 = FiberSingle(12, -4, 14, 2, 1)
a3 = FiberSingle(12, 4, -14, 2, 1)
a4 = FiberSingle(12, 4, 14, 2, 1)
c.add_fibers(b1, a1, a2, a3, a4)

# Output Information From FiberSection
c.print_section_properties()
c.plot_fibers()

# Build ACI Strain Compatibility Object
test = AciStrainCompatibility(c)
test.add_concrete_boundary(0.5 * B, 0.5 * H, 0)
test.add_concrete_boundary(0.5 * B, -0.5 * H, 0)
test.add_concrete_boundary(-0.5 * B, -0.5 * H, 0)
test.add_concrete_boundary(-0.5 * B, 0.5 * H, 0)
test.add_steel_boundary(4, -14, 0)
test.add_steel_boundary(4, 14, 0)
test.add_steel_boundary(-4, 14, 0)
test.add_steel_boundary(-4, -14, 0)
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
