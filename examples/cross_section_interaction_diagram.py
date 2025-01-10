from libdenavit import CrossSection2d
from libdenavit.section import RC, Circle, ReinfCirc
from math import pi
import matplotlib.pyplot as plt

# Input Properties
length = 240
D = 16
num_bars = 8
fc = 4
fy = 60
rhosr = 0.02
Ab = rhosr * pi / 4 * D ** 2 / num_bars
axis = 'x'
start_material_id = 1

# Define RC section object
conc_cross_section = Circle(D)
reinforcement = ReinfCirc(0.5 * D - 1, num_bars, Ab)
section = RC(conc_cross_section, reinforcement, fc, fy, 'US')
# section.plot_section()

col = CrossSection2d(section, axis=axis)
section_args = (start_material_id, "ElasticPP", "Concrete04_no_confinement", 20, 20)
section_kwargs = dict()

result_ops = col.run_ops_interaction(section_args=section_args, section_kwargs=section_kwargs, 
                                     prop_load_incr_factor=1e-2, nonprop_load_incr_factor=1e-1)
result_design = col.run_AASHTO_interaction(section_factored=False)
plt.plot(result_ops["M1"], result_ops["P"], '-', label='OPS')
plt.plot(result_design["M1"], -result_design["P"], 'r-', label='Design')
plt.xlabel('Bending Moment (kip-in.)')
plt.ylabel('Axial Force (kips)')
plt.xlim(0, )
plt.ylim(0, )
plt.legend()
plt.show()
