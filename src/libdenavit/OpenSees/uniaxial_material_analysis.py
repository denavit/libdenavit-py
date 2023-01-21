import openseespy.opensees as ops
import matplotlib.pyplot as plt
from math import ceil

def uniaxial_material_analysis(definition, peak_points, rate_type = 'None', rate_value = None, 
                               def_args = [], def_kwargs = {}, matTag=1, ndm = 1, ndf = 1, parallel_stiffness = 0,
                               plot_stress_strain = False, compression_positive = False):
        
    # ------------------------------
    # Build Model
    # -----------------------------        
    ops.wipe()
    ops.model('basic', '-ndm', ndm, '-ndf', ndf)
    
    # Define Nodes
    if ndm == 1:
        ops.node(1,0.0)
        ops.node(2,1.0)
    elif ndm == 2:
        ops.node(1,0.0,0.0)
        ops.node(2,1.0,0.0)
    elif ndm == 3:
        ops.node(1,0.0,0.0,0.0)
        ops.node(2,1.0,0.0,0.0)
    else:
        raise ValueError('ndm should be 1, 2, or 3')
    
    # Define Boundary Conditions
    if ndf == 1:
        ops.fix(1, 1)
    elif ndf == 3:
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 1)
    else:
        raise ValueError(f'Code not implemented for ndf = {ndf}')

    # Define Material
    if isinstance(definition, list):
        ops.uniaxialMaterial(*definition)
    elif callable(definition):
        definition(*def_args,**def_kwargs)
    else:
        raise TypeError('Material defintion is not a type that can be used')
       
    # Define Element
    ops.element('Truss', 1, 1, 2, 1.0, matTag)
    
    # Add parallel stiffness (if requested)
    if parallel_stiffness > 0:
        ops.uniaxialMaterial('Elastic', 24325, parallel_stiffness)
        ops.element('Truss', 2, 1, 2, 1.0, 24325)
        
    # Define Loads
    ops.timeSeries('Linear', 854)
    ops.pattern('Plain', 979, 854)
    if ndf == 1:
        ops.load(2, 1.0)
    elif ndf == 3:
        ops.load(2, 1.0, 0.0, 0.0)
    else:
        raise ValueError(f'Code not implemented for ndf = {ndf}')

    # ------------------------------
    # Build and Perform the Analysis
    # ------------------------------       
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormUnbalance', 1e-8, 10)
    ops.algorithm('Newton')
    ops.integrator('DisplacementControl', 2, 1, peak_points[0])
    ops.analysis('Static')
    ops.analyze(1)
    if compression_positive:
        strain = [-1*ops.nodeDisp(2,1)]
        stress = [-1*ops.getTime() - parallel_stiffness*strain[-1]]
    else:
        strain = [ops.nodeDisp(2,1)]
        stress = [ops.getTime() + parallel_stiffness*strain[-1]]
    
    for i in range(len(peak_points)-1):
        if rate_type == 'None':
            num_steps = 1
        elif rate_type == 'StrainRate':
            num_steps = ceil(abs(peak_points[i+1] - peak_points[i])/rate_value)
        elif rate_type == 'Steps':
            num_steps = rate_value
        else:
            raise ValueError('Unknown rate_type: %s' % rate_type)
            
        incr = (peak_points[i+1] - peak_points[i])/num_steps
        
        for j in range(num_steps):
            ops.integrator('DisplacementControl', 2, 1, incr)
            ops.analyze(1)
            if compression_positive:
                strain.append(-1*ops.nodeDisp(2,1))
                stress.append(-1*ops.getTime() - parallel_stiffness*strain[-1])
            else:
                strain.append(ops.nodeDisp(2,1))
                stress.append(ops.getTime() + parallel_stiffness*strain[-1])
            
    if plot_stress_strain:
        fig = plt.figure()
        plt.plot(strain,stress)
        plt.xlabel('Strain')
        plt.ylabel('Stress')
        plt.show()
            
    return (strain,stress)