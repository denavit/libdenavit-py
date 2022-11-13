import openseespy.opensees as ops
import matplotlib.pyplot as plt
from math import ceil

def uniaxial_material_analysis(definition, peak_points, rate_type = 'None', rate_value = None, matTag=1, 
                               plot_stress_strain = False, compression_positive = False):
        
    # ------------------------------
    # Build Model
    # -----------------------------        
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Define Nodes
    ops.node(1,0.0)
    ops.node(2,1.0)

    # Define Boundary Conditions
    ops.fix(1, 1)

    # Define Elements
    ops.uniaxialMaterial(*definition)
    ops.element('Truss', 1, 1, 2, 1.0, matTag)
     
    # Define Loads
    ops.timeSeries('Linear', 854)
    ops.pattern('Plain', 979, 854)
    ops.load(2, 1.0)

    # ------------------------------
    # Build and Perform the Analysis
    # ------------------------------       
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandSPD')
    ops.test('NormUnbalance', 1e-8, 10)
    ops.algorithm('Newton')
    ops.integrator('DisplacementControl', 2, 1, peak_points[0])
    ops.analysis('Static')
    ops.analyze(1)
    if compression_positive:
        stress = [-1*ops.getTime()]
        strain = [-1*ops.nodeDisp(2,1)]
    else:
        stress = [ops.getTime()]
        strain = [ops.nodeDisp(2,1)]
    
    for i in range(len(peak_points)-1):
        if rate_type == 'None':
            num_steps = 1
        elif rate_type == 'StrainRate':
            num_steps = ceil(abs(peak_points[i+1] - peak_points[i])/rate_value)
        elif rate_type == 'Steps':
            num_steps = rate_value
        else:
            raise Exception('Unknown rate_type: %s' % rate_type)
            
        incr = (peak_points[i+1] - peak_points[i])/num_steps
        
        for j in range(num_steps):
            ops.integrator('DisplacementControl', 2, 1, incr)
            ops.analyze(1)
            if compression_positive:
                stress.append(-1*ops.getTime())
                strain.append(-1*ops.nodeDisp(2,1))
            else:
                stress.append(ops.getTime())
                strain.append(ops.nodeDisp(2,1))
            
    if plot_stress_strain:
        fig = plt.figure()
        plt.plot(strain,stress)
        plt.xlabel('Strain')
        plt.ylabel('Stress')
        plt.show()
            
    return (strain,stress)