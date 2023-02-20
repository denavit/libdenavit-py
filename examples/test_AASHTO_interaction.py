from libdenavit import NonSwayColumn2d, InteractionDiagram2d
from libdenavit.section import RC, Circle, ReinfCirc
import matplotlib.pyplot as plt
import math
import numpy as np

# Input Properties
D = 48
num_bars = 24
fc = 4
fy = 60
rhosr = 0.02
Ab = rhosr * math.pi/4*D**2 / (num_bars)
length = 15 * D
axis = 'x'

# Define RC section object
conc_cross_section = Circle(D)
reinforcement = ReinfCirc(0.5*D - 3, num_bars, Ab)
section = RC(conc_cross_section, reinforcement, fc, fy, 'US')

# Plot cross section
section.plot_section(show=False)

# Plot cross-sectional interaction diagram 
P_CS1, M_CS1, _ = section.section_interaction_2d(axis, 100, factored=False)
P_CS2, M_CS2, _ = section.section_interaction_2d(axis, 100, factored=True)
fig = plt.figure()
plt.plot(M_CS1, -P_CS1, color='tab:orange', linestyle='--')
plt.plot(M_CS2, -P_CS2, color='tab:orange', linestyle='--')
plt.xlabel("Bending Moment (kip-in.)")
plt.ylabel("Axial Compression (kips)")

# Create NonSwayColumn2d object
col = NonSwayColumn2d(section, length, 1, 1, dxo=length / 1000, axis=axis)

# Define fiber section definition arguments 
section_args = (1, "ElasticPP", "Concrete04_no_confinement", 20, 20)
section_kwargs = {}

# Run Analyses
num_points = 10
Opensees_results = col.run_ops_interaction(section_args, section_kwargs, num_points)
AASHTO_results  = col.run_AASHTO_interaction(EI_type='a', section_factored = False, Pc_factor = 1.00)
# Make Plot
fig = plt.figure()
plt.plot(Opensees_results["M1"], Opensees_results["P"], color='tab:blue',   marker='x', linestyle='--', label='$M_1$ (OpenSees)')
plt.plot(Opensees_results["M2"], Opensees_results["P"], color='tab:blue',   marker='o', linestyle='-',  label='$M_2$ (OpenSees)')
plt.plot(AASHTO_results["M1"], AASHTO_results["P"],  color='tab:orange', marker='x', linestyle='--', label='$M_1$ (AASHTO)')
plt.plot(AASHTO_results["M2"], AASHTO_results["P"],  color='tab:orange', marker='o', linestyle='-',  label='$M_2$ (AASHTO)')
plt.xlabel("Bending Moment (kip-in.)")
plt.ylabel("Axial Compression (kips)")
plt.legend(loc='upper right')


angles = [x + 0.1 for x in range(0, 91)]
angles[-1] = angles[-1] - 0.2

AASHTO = InteractionDiagram2d(AASHTO_results["M1"], AASHTO_results["P"], is_closed=True)
Opensees = InteractionDiagram2d(Opensees_results["M1"], Opensees_results["P"], is_closed=True)
e = Opensees.compare_two(AASHTO, angles, degrees=True)

fig = plt.figure(figsize=(6, 4))
plt.plot(angles, e)
plt.axhline(y=0, color='k')
plt.ylabel(r'$\epsilon= \frac{r_{OpenSees} - r_{AASHTO}}{r_{OpenSees}}$')
plt.xlabel('Angle with respect to x-axis (degrees)')
plt.tight_layout()
plt.xlim(0, 90)
plt.xticks(np.arange(0, 91, 10))

# Show Plots
plt.show()
