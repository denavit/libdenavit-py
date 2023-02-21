from libdenavit import SwayColumn2d, NonSwayColumn2d
from libdenavit.section import RC, Rectangle, ReinfRect
from math import pi, inf
import matplotlib.pyplot as plt
import numpy as np

result_ops = None
result_design = None
nonsway = False
sway = True
plot_section = False

# region Input Properties
length = 11.25 * 12
fc = 4
fy = 60
rhosr = 0.02
Ab = 0.3
axis = 'x'
G_top = inf
G_bot = 1.48
# endregion

# region Define RC section object
conc_cross_section = Rectangle(12, 12)
reinforcement = ReinfRect(10, 10, 2, 2, Ab)
section = RC(conc_cross_section, reinforcement, fc, fy, 'US')
if plot_section:
    section.plot_section()
EIeff = section.EIeff(axis, "b", 0)
Ec = section.Ec
print(f"{Ec= }")
print(f"{EIeff= }")
# endregion

# region Non-sway
if nonsway:
    col = NonSwayColumn2d(section, length, 1, 1, axis=axis)
# endregion

# region Sway definition
elif sway:
    Igc = section.EIgross(axis)
    k_top = 6 * (0.4 * Ec * Igc) / (G_top * length)
    k_bot = 6 * (0.4 * Ec * Igc) / (G_bot * length)
    col = SwayColumn2d(section, length, k_bot, k_top, 0, axis=axis)
    K = col.effective_length_factor(EIeff)
    print(f"{K= }")

    col._K = K

    Pc = pi ** 2 * EIeff / (K * col.length) ** 2
    print(f'{Pc= }')
    print(f'{col.Cm= }')
    P = 66
    delta_s = 1 / (1 - P / (0.75 * Pc))
    print(f'{delta_s= }')
# endregion

# region run OPS
# section_args = (1, "ElasticPP", "Concrete04_no_confinement", 20, 20)
# section_kwargs = dict()
#
# result = col.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=0)
# results = col.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=0,
#                                         disp_incr_factor=1e-4)
# print(f"P: {result.applied_axial_load}")
# print(f"M_top: {result.moment_at_top}")
# print(f"M_bottom: {result.moment_at_bottom}")
# print(f"Max abs disp: {result.maximum_abs_disp}")
# print(result.exit_message)
# endregion

# region run interaction
# result_ops = col.run_ops_interaction(section_args, section_kwargs)
result_design = col.run_AASHTO_interaction('b', section_factored=False)
# endregion

# region plot and print results
if result_ops is not None:
    print(f"{result_ops['M1']= }")
    print(f"{result_ops['P']= }")
    plt.plot(result_ops["M1"], result_ops["P"], 'b-', label='ops')
if result_design is not None:
    plt.plot(result_design["M1"], np.array(result_design["P"]), 'r-', label='Design')
P_CS, M_CS, _ = col.section.section_interaction_2d(col.axis, 100, factored=False, only_compressive=True)
plt.plot(M_CS, P_CS, 'k-', label='Section')
plt.legend()
plt.show()
# endregion

# region calculate EI
print(f"{result_design['P']= }")
results_EI = col.calculated_EI(result_design["P"], result_design["M1"], P_CS=P_CS, M_CS=M_CS)
print(f'P: {results_EI["P"]}')
print(f'Calculated EI: {results_EI["EI"]}')
print(f'EI gross: {results_EI["EIgross"]}')
# endregion

