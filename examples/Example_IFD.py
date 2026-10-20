import openseespy.opensees as ops
import libdenavit.section.database.aisc as section
from libdenavit.OpenSees.plotting import plot_undeformed_2d,plot_deformed_2d,plot_element_values_2d,plot_internal_force_diagram_2d

E=29000
L_x=5*12
H=3*12
Node_Details=[[1,0,0],[2,0,H],[3,L_x,H],[4,L_x,0]]
Connectivity_Details= [[1,2,'W8X31'], [2,3,'W8X31'],[3,4,'W8X31']]  ## [[nodei1,nodej1,section1],[]]


ops.wipe()
ops.model('basic','-ndm',2,'-ndf',3)
for single_node in Node_Details:
    ops.node(single_node[0],single_node[1],single_node[2])

ops.fix(1,1,1,1)
ops.fix(4,1,1,1)
ops.geomTransf("Linear", 1)

# Define Elastic Elements with Connectivity Details
for i,element_details in enumerate(Connectivity_Details):
    element_node1=element_details[0]
    element_node2=element_details[1]
    section_name=element_details[2]
    A=section.wide_flange_database[section_name]['A']
    Ix=section.wide_flange_database[section_name]['Ix']
    ops.element('elasticBeamColumn',i+1,element_node1,element_node2,A,E,Ix,1)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 50,0,0)
ops.load(3,0,-30,0)
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('UmfPack')
ops.integrator('LoadControl', 1/10)
ops.algorithm('Linear')
ops.analysis('Static')
ops.analyze(10)

plot_undeformed_2d()
plot_deformed_2d(scale_factor=50)
plot_internal_force_diagram_2d(force_type='P',show_values=True)
plot_internal_force_diagram_2d(force_type='V',scale=0.3,show_values=True)
plot_internal_force_diagram_2d(force_type='M',scale=0.02,show_values=True)













