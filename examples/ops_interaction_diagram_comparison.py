from libdenavit import NonSwayColumn2d
from libdenavit.section import RC, Rectangle, ReinfRect, Circle
import matplotlib.pyplot as plt

# Input Properties
H = 20
B = 15
nbH = 2
nbB = 2
lh = 10
fc = 4
fy = 60
rhosr = 0.02
Ab = rhosr * H * B / (nbH * nbB - 2)
length = lh * H
axis = 'x'

# Define RC section object
conc_cross_section = Rectangle(H, B)
reinforcement = ReinfRect(H - 0.15 * H, B - 0.15 * H, nbB, nbH, Ab)
section = RC(conc_cross_section, reinforcement, fc, fy, 'US')

# Create NonSwayColumn2d object
col = NonSwayColumn2d(section, length, 1, 1, dxo=length / 1000)

# Define fiber section definition arguments 
section_args = (1, "ElasticPP", "Concrete04_no_confinement", 20, 20)
section_kwargs = {'axis':axis}

# Run Analyses
num_points = 10
P_np, M1_np, M2_np = col.run_ops_interaction(section_args, section_kwargs, num_points)

h = section.depth(axis) 
e_list = [0,0.01*h,0.02*h,0.04*h,0.06*h,0.08*h,0.1*h,0.2*h,0.4*h,0.7*h,1.0*h]
P_p,  M1_p,  M2_p  = col.run_ops_interaction_proportional(section_args, section_kwargs, e_list, disp_incr_factor=0.00001)

# Make Plot
fig = plt.figure()
plt.plot(M1_np, P_np, color='tab:blue',   marker='x', linestyle='--', label='$M_1$ (non-prop.)')
plt.plot(M2_np, P_np, color='tab:blue',   marker='o', linestyle='-',  label='$M_2$ (non-prop.)')
plt.plot(M1_p,  P_p,  color='tab:orange', marker='x', linestyle='--', label='$M_1$ (prop.)')
plt.plot(M2_p,  P_p,  color='tab:orange', marker='o', linestyle='-',  label='$M_2$ (prop.)')
plt.xlabel("Bending Moment (kip-in.)")
plt.ylabel("Axial Compression (kips)")
plt.legend(loc='upper right')
plt.show()
