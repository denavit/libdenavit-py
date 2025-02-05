import openseespy.opensees as ops
import matplotlib.pyplot as plt
import math

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

N = 1.0
mm = 1.0

MPa = N / mm ** 2
GPa = 1000 * MPa
kN = 1000 * N

b = 300 * mm
h = b

Ag = b * h
As = 1500 * mm ** 2
Ac = Ag - As

P = 1000 * kN
Ec = 25 * GPa
Es = 200 * GPa

fc = ((Ec / MPa) / 4700) ** 2 * MPa
epsshu0 = 600e-6
phiu0 = 3.0

A = b * h
O = 2 * b + 2 * h
VoverS = A / O
h0 = 2 * VoverS
psish = 26 * math.exp(0.0142 * VoverS)
phiu = 3.0
psicr1 = 1.0
psicr2 = 75.4
tcast = 0
tDry = 14
Tcr = 28


ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(1, 0, 0)
ops.fix(1, 1, 1, 1)
ops.node(2, 0, 0)
ops.fix(2, 0, 1, 1)

# Creep wrapper with elastic concrete
ops.uniaxialMaterial('Elastic', 2, Ec)
ops.uniaxialMaterial('Creep', 22, 2, tDry, -epsshu0, psish, Tcr, phiu, psicr1, psicr2, tcast)
# Creep wrapper with Concrete02
ops.uniaxialMaterial('Concrete02', 3, fc, 0.002, 0, 0.006)
ops.uniaxialMaterial('Creep', 33, 3, tDry, -epsshu0, psish, Tcr, phiu, psicr1, psicr2, tcast)

concTag = 22

ops.section('Fiber', 1)
ops.fiber(0, 0, Ac, concTag)

ops.element('zeroLengthSection', 1, 1, 2, 1)

# Static analysis of zero load
ops.setTime(tDry)
#ops.test('NormDispIncr',1e-3,10,0)
ops.integrator('LoadControl', 0.0)
ops.setCreep(0)
ops.system('UmfPack')
ops.analysis('Static', '-noWarnings')
ops.analyze(1)

# Gravity load
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, -P, 0, 0)

# Static analysis of gravity load
ops.setTime(tDry)
#ops.test('NormDispIncr',1e-3,10,0)
ops.integrator('LoadControl', 0.0)
ops.setCreep(0)
ops.analysis('Static', '-noWarnings')
ops.analyze(1)

# Stresses at start of creep
sigc = ops.eleResponse(1, 'section', 'fiber', 0, 0, concTag, 'stress')[0] / MPa
eps_c = ops.eleResponse(1, 'section', 'fiber', 0, 0, concTag, 'strain')[0]
cplot = [-sigc]

eps_c_plot = [eps_c]

ops.setTime(tDry)
ops.setCreep(1)

t = Tcr

tplot = [t]

# Start and end times (days)
t0 = 28
tfinish = 1000000

t = t0
while t < tfinish:
    logt0 = math.log10(t0)
    logt1 = logt0 + 0.001
    t1 = 10 ** logt1
    dt = t1 - t0  # Change in time (days)
    t0 = t1
    ops.integrator('LoadControl', dt)
    t = t1

    ok = ops.analyze(1)
    if ok < 0:
        break

    sigc = ops.eleResponse(1, 'section', 'fiber', 0, 0, concTag, 'stress')[0] / MPa

    epsc = ops.eleResponse(1, 'section', 'fiber', 0, 0, concTag, 'strain')[0]

    cplot.append(-sigc)
    tplot.append(t)

    eps_c_plot.append(epsc)

fig2, ax3 = plt.subplots()
ax3.semilogx(tplot, eps_c_plot, color='tab:blue', label='Concrete')
ax3.tick_params(axis='y',)
ax3.grid()
ax3.set_ylabel('Concrete Strain')
ax3.set_xlabel('Time (days)')
ax3.set_xlim(right=tfinish)
plt.tight_layout()
plt.show()
