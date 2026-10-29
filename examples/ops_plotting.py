from libdenavit import opensees as ops
from libdenavit.section.database.aisc import wide_flange_database
from libdenavit.OpenSees.plotting import plot_undeformed_2d,plot_deformed_2d,plot_element_values_2d,plot_internal_force_diagram_2d, animate_element_values_2d
import matplotlib.pyplot as plt

inch = 1
ft = 12*inch
kip = 1
ksi = kip/inch**2


ops.wipe()
ops.model('basic','-ndm',2,'-ndf',3)

# Define nodes
ops.node(1,  0.0,  0.0)
ops.node(2,  0.0,12*ft)
ops.node(3,20*ft,12*ft)
ops.node(4,20*ft,  0.0)

ops.fix(1,1,1,1)
ops.fix(4,1,1,1)

# Define elements
E = 29000*ksi
column_section = wide_flange_database['W8X31']
beam_section = wide_flange_database['W18X35']

ops.geomTransf("Linear", 1)
ops.element('elasticBeamColumn',1,1,2,column_section['A'],E,column_section['Ix'],1)
ops.element('elasticBeamColumn',2,2,3,beam_section['A'],E,beam_section['Ix'],1)
ops.element('elasticBeamColumn',3,3,4,column_section['A'],E,column_section['Ix'],1)

# Define loads
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 10.0*kip,-30*kip,0)
ops.load(3,        0,-30*kip,0)

# Define and run analysis
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('UmfPack')
ops.integrator('LoadControl', 1)
ops.algorithm('Linear')
ops.analysis('Static')
ops.analyze(1)


# Make plots
plot_undeformed_2d()

plot_deformed_2d(scale_factor=50)

element_data = {1:0.10, 2:0.30, 3:0.50}
plot_element_values_2d(element_data)

plot_internal_force_diagram_2d(force_type='P',scale=0.6,annotate_max_min=False,annotate=True)

plot_internal_force_diagram_2d(force_type='V',scale=1,annotate_max_min=False,annotate=True)

plot_internal_force_diagram_2d(force_type='M',scale=0.06,annotate_max_min=False,annotate=True,new_figure=False)

element_data_history = [
    {1: 0.10, 2: 0.30, 3: 0.50},
    {1: 0.15, 2: 0.35, 3: 0.55},
    {1: 0.20, 2: 0.40, 3: 0.60},
    {1: 0.25, 2: 0.45, 3: 0.65},
    {1: 0.30, 2: 0.50, 3: 0.70},
    {1: 0.35, 2: 0.55, 3: 0.75},
    {1: 0.40, 2: 0.60, 3: 0.80},
    {1: 0.45, 2: 0.65, 3: 0.85}
]

animate_element_values_2d(element_data_history,interval=200)













