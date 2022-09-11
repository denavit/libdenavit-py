from libdenavit import NonSwayColumn2d
from libdenavit.section import RC, Circle, ReinfCirc
import matplotlib.pyplot as plt
import math

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
col = NonSwayColumn2d(section, length, 1, 1, dxo=length / 1000)

# Define fiber section definition arguments 
section_args = (1, "ElasticPP", "Concrete04_no_confinement", 20, 20)
section_kwargs = {'axis':axis}

# Run Analyses
num_points = 10
P_ops, M1_ops, M2_ops = col.run_ops_interaction(section_args, section_kwargs, num_points)

P_AASHTO, M1_AASHTO, M2_AASHTO  = col.run_AASHTO_interaction(axis, EI_type='a', section_factored = False, Pc_factor = 1.00)


# Make Plot
fig = plt.figure()
plt.plot(M1_ops, P_ops, color='tab:blue',   marker='x', linestyle='--', label='$M_1$ (OpenSees)')
plt.plot(M2_ops, P_ops, color='tab:blue',   marker='o', linestyle='-',  label='$M_2$ (OpenSees)')
plt.plot(M_CS1, -P_CS1, color='tab:orange', linestyle='-.')
plt.plot(M1_AASHTO, -P_AASHTO,  color='tab:orange', marker='x', linestyle='--', label='$M_1$ (AASHTO)')
plt.plot(M2_AASHTO, -P_AASHTO,  color='tab:orange', marker='o', linestyle='-',  label='$M_2$ (AASHTO)')
plt.xlabel("Bending Moment (kip-in.)")
plt.ylabel("Axial Compression (kips)")
plt.legend(loc='upper right')


# Show Plots
plt.show()
