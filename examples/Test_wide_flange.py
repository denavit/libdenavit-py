import math as math
from libdenavit.cross_section_2d import *
from libdenavit.section.wide_flange import *
import matplotlib.pyplot as plt


#################################################
# E=29000*ksi
# fy=36*ksi
# Hk = 0.001*E           # Kinematic hardening modulus
E=199947961
fy=248210
Hk = 0.001*E

#Values for W27X84
d=0.678
tw=0.011
bf=0.254
tf=0.016
A=0.016
Ix=0.001186
Zx=0.00399
Iy=4.41205*10**(-5)
Zy=0.000544

beam_section_tag=1
axis='x'
Mp=fy*Zx
Initial_slope=E*Ix

def Moment_curvature_analysis(P_value, residual_stress, axis):
    frc = -0.3 * fy if residual_stress else 0
    beam = I_shape(d, tw, bf, tf,
                   Fy=fy, E=E, Hk=Hk,
                   A=A, Ix=Ix, Iy=Iy)
    member = CrossSection2d(beam, axis=axis)
    results = member.run_ops_analysis(
        analysis_type='nonproportional_limit_point',
        **{
            'section_id': beam_section_tag,
            'section_args': [1, 'Steel01', 20, 20, frc],
            'section_kwargs': {},
            'P': P_value
        }
    )
    return results

def axial_load_comparison(P_values=(0, 500, 1000, 1500, 2000), residual_stress=True, axis=axis):
    plt.figure()
    for P in P_values:
        results = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=axis)
        label = f"P={results.applied_axial_load[-1]:.0f}"
        curv = results.curvatureX if axis == 'x' else results.curvatureY
        plt.plot(curv, results.maximum_abs_moment, label=label)

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature vs Axial Load ({rs_text})")
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_with_EI_and_FyZ(P=0, residual_stress=False, axis='x'):
    plt.figure()
    results = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=axis)
    label = f"P={results.applied_axial_load[-1]:.0f}"
    curv = results.curvatureX if axis == 'x' else results.curvatureY
    plt.plot(curv, results.maximum_abs_moment, label=label)

    plt.axhline(y=Mp, linestyle='--', linewidth=1.2, label=f"M_p = {Mp:.3g}")
    kappa_end = Mp / Initial_slope
    plt.plot([0, kappa_end], [0, Mp], linestyle=':', linewidth=1.2, label=f"slope = {Initial_slope:.3g}")

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature vs Axial Load ({rs_text})")
    plt.grid(True)
    plt.legend()
    plt.show()

def residual_stress_comparison():
    P_values = [0,  1000]
    rs_flags = [False, True]  # without RS, with RS
    styles = {False: {'linestyle': '-',  'linewidth': 1.4},
            True:  {'linestyle': '--', 'linewidth': 1.4}}
    plt.figure()

    for P in P_values:
        for rs in rs_flags:
            results = Moment_curvature_analysis(P_value=P, residual_stress=rs, axis=axis)
            label = f"P={results.applied_axial_load[-1]:.0f}, RS={'Yes' if rs else 'No'}"
            curv = results.curvatureX if axis == 'x' else results.curvatureY
            plt.plot(curv, results.maximum_abs_moment, label=label, **styles[rs])

    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title("Moment–Curvature Curves (with/without Residual Stress)")
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_of_major_and_minor_axes(P=0, axes=('x', 'y')):
    plt.figure()
    for rs in [False, True]:  # Loop over residual stress states
        rs_label = "with RS" if rs else "without RS"
        linestyle = '-' if not rs else '--'  # solid for no RS, dashed for with RS

        for ax in axes:
            res = Moment_curvature_analysis(P_value=P, residual_stress=rs, axis=ax)
            curv = res.curvatureX if ax == 'x' else res.curvatureY

            plt.plot(curv, res.maximum_abs_moment,
                     linestyle=linestyle,
                     label=f"{ax}-axis ({rs_label})")

            # Elastic-plastic line
            Mp_axis = fy * (Zx if ax == 'x' else Zy)
            EI_axis = E * (Ix if ax == 'x' else Iy)
            kappa_end = Mp_axis / EI_axis
            if not rs:  # Only plot the slope line once (for uncluttered view)
                plt.axhline(y=Mp_axis, linestyle=':', linewidth=1.2,
                            label=f"M_p,{ax} = {Mp_axis:.3g}")
                plt.plot([0, kappa_end], [0, Mp_axis], linestyle=':', linewidth=1.2,
                         label=f"slope EI_{ax} = {EI_axis:.3g}")

    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature Comparison (P = {P})")
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_of_2d_and_3d_section(P=0, residual_stress=False, axes=('x','y',None)):
    plt.figure()
    for ax in axes:
        res = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=ax)
        curv = res.curvatureX if ax == 'x' or ax== None else res.curvatureY
        plt.plot(curv, res.maximum_abs_moment, label=f"P={res.applied_axial_load[-1]:.0f}, axis={ax}")

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature (axes={list(axes)}, {rs_text}, P={P})")
    plt.grid(True)
    plt.legend()
    plt.show()

residual_stress_comparison()
axial_load_comparison()
comparison_with_EI_and_FyZ()
comparison_of_major_and_minor_axes()
comparison_of_2d_and_3d_section()