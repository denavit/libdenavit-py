from math import sqrt
from libdenavit.OpenSees import get_fiber_data
import openseespy.opensees as ops

'''
This is an example showing how to use libdenavit.OpenSees.get_fiber_data
'''

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.node(0, 0, 0, 0)
ops.node(1, 0, 0, 0)

ops.fix(1, 1, 0, 1, 1, 0, 0)
ops.uniaxialMaterial('Elastic', 1, 1000)

# Define steel angle cross section
b_angle = 3
t_angle = 1/4

ops.section('Fiber', 1, '-GJ', 10e6)
yI = -b_angle/sqrt(2) - t_angle/sqrt(2)
zI = -b_angle/sqrt(2) + t_angle/sqrt(2)
yJ = -b_angle/sqrt(2)
zJ = -b_angle/sqrt(2)
yK = 0.0
zK = 0.0
yL = -sqrt(2)*t_angle
zL = 0.0
ops.patch('quad', 1, 3, 20, yI, zI, yJ, zJ, yK, zK, yL, zL)
yI = -sqrt(2)*t_angle
zI = 0.0
yJ = 0.0
zJ = 0.0 
yK = -b_angle/sqrt(2)
zK =  b_angle/sqrt(2)
yL = -b_angle/sqrt(2) - t_angle/sqrt(2)
zL =  b_angle/sqrt(2) - t_angle/sqrt(2)
ops.patch('quad', 1, 3, 20, yI, zI, yJ, zJ, yK, zK, yL, zL)


ops.element('zeroLengthSection', 1, 0, 1, 1)

x, y, A, m = get_fiber_data(1, plot_fibers = True)
print(f'x: {x}\n', f'y: {y}\n', f'A: {A}\n', f'm: {m}\n')