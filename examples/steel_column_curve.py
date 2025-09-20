from libdenavit.section import I_shape
from libdenavit import NonSwayColumn2d
import numpy as np
import matplotlib.pyplot as plt



def run_single_analysis(L,residual_stress,material_options):
    """
    Run nonlinear column analysis for a given length L
    Returns: Peak axial load, moment, and displacement at limit point
    """

    section_kwargs = dict(
        start_material_id=1,
        nfy=20,
        nfx=20,
        frc=residual_stress,
        stiffness_reduction=0.9,
        strength_reduction=0.9,
        GJ=1.0e6,
    )
    section_kwargs.update(material_options)

    col = NonSwayColumn2d(
        section=sec,
        length=L,
        et=0.,
        eb=0.,
        axis='y',
        dxo=L/1000.0,
        n_elem=8,
        element_type='mixedBeamColumn',
        ops_geom_transf_type='Corotational',
        ops_integration_points=3,
    )

    analysis_kwargs = dict(
        section_id=1,
        section_args=[],
        section_kwargs=section_kwargs,
        e=0.0,
        disp_incr_factor=1e-5,
        percent_load_drop_limit=0.05,
        steel_strain_limit=0.05,
        deformation_limit=0.05*L,
    )    

    res = col.run_ops_analysis("proportional_limit_point", **analysis_kwargs)

    Pstar = getattr(res, "applied_axial_load_at_limit_point", np.nan)
    Mstar = getattr(res, "maximum_abs_moment_at_limit_point", np.nan)
    Dstar = getattr(res, "maximum_abs_disp_at_limit_point", np.nan)

    print(f"L = {L:.1f} in | P* = {Pstar:.6f} kips | M* = {Mstar:.6f} kip-in | Δ* = {Dstar:.6f} in | exit = {getattr(res, 'exit_message', 'None')}")
    return Pstar, Mstar, Dstar




# === Input Data ===
section_name = 'W14x145'
E  = 29000
Fy = 50

sec = I_shape.from_database(section_name, Fy=Fy, E=E)
A = sec.A
Py = A*Fy

length_range = np.linspace(50,500,19)

Steel01_options = dict(
    mat_type="Steel01",
    hardening_ratio=0.001,
)

MultiLinear_options = dict(
    mat_type='MultiLinear',
)

analysis_options_name_list = ['Steel01 with RS','Steel01 without RS','MultiLinear']
residual_stress_list = [-0.3*Fy,0,0]
material_options_list = [Steel01_options,Steel01_options,MultiLinear_options]



# === MAIN ANALYSIS LOOP ===
plt.figure()

for analysis_options_name,residual_stress,material_options in zip(analysis_options_name_list,residual_stress_list,material_options_list):
    P_cap, M_cap, D_cap = [], [], []

    for L in length_range:
        Pstar, Mstar, Dstar = run_single_analysis(L,residual_stress,material_options)
        P_cap.append(Pstar)
        M_cap.append(Mstar)
        D_cap.append(Dstar)

    plt.plot(length_range, np.array(P_cap)/Py, 'o-', lw=2, label=analysis_options_name)

plt.legend(loc='upper right')
plt.xlim(0,max(length_range))
plt.ylim(0,1.1)
plt.xlabel("Column Length, L (in.)")
plt.ylabel("Normalized Axial Strength, $P/A_gF_y$")
plt.grid(True)
plt.show()
