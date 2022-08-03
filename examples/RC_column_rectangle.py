import matplotlib.pyplot as plt
from libdenavit.section import FiberSection, FiberQuadPatch, FiberSingle, AciStrainCompatibility
import numpy as np


# A simple example of making a fiber section (This one does not use RC class)

# Units: kips, in, ksi
B  = 20     # Cross-sectional width
H  = 40     # Cross-sectional height
fc = 4      # Concrete compressive strength
fy = 60     # Steel yield strength
Es = 29000
cover = 0.15 * H
rhosr = 0.06
nbB = 5
nbH = 2
axis = 'x'
units = 'us'
Ab = H*B*rhosr/(2*nbB+2*nbH-4)

# Build Fiber Section Object
c = FiberSection(100, 100)
b1 = FiberQuadPatch(-10, -20, -10, 20, 10, 20, 10, -20, 1)  # xI, yI, ..., m
a1 = FiberSingle(Ab, -4, -14, 2, 1)  # A, x, y, m, m_neg
a2 = FiberSingle(Ab, -2, -14, 2, 1)
a4 = FiberSingle(Ab, -0, -14, 2, 1)
a5 = FiberSingle(Ab,  2, -14, 2, 1)
a3 = FiberSingle(Ab,  4, -14, 2, 1)
a6  = FiberSingle(Ab, -4, 14, 2, 1)
a7  = FiberSingle(Ab, -2, 14, 2, 1)
a8  = FiberSingle(Ab,  0, 14, 2, 1)
a9  = FiberSingle(Ab,  2, 14, 2, 1)
a10 = FiberSingle(Ab,  4, 14, 2, 1)

c.add_fibers(b1, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)

# Output Information From FiberSection
properties = c.print_section_properties()
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
b = np.arange(-30, 20, 0.5)
c = np.arange(20, 80, 0.5)
d = np.concatenate((a, b, c))

# Perform Interaction Calculations
plot_P = []
plot_Mx = []
for i in d:
    P, Mx, My, et = test.compute_point(0, i, 0)
    plot_P.append(-P)
    plot_Mx.append(Mx)

Ag = properties['Area'].iloc[-1]

plot_Mx = [x/ (fc * B * H * H) for x in plot_Mx]
plot_P = [x / (fc * B * H) for x in plot_P]

# Plot Interaction Diagram
plt.xlim(0, 0.55)
plt.ylim(0, 2)
plt.xticks(np.arange(0, 0.60, 0.05))
plt.ylabel(r"$K_n=P_n/f'_cA_g$")
plt.xlabel(r"$K_x=P_ne/f'_cA_gh$")
plt.title("Interaction Diagram")
plt.plot(plot_Mx, plot_P, 'b-')
plt.show()
