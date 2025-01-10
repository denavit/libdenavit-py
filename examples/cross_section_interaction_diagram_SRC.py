from libdenavit import CrossSection2d
from libdenavit.section import RC, Circle, ReinfCirc
from math import pi
import matplotlib.pyplot as plt

# Input Properties

shape = 'W14X99'

Fy = 50
Fylr = 60
fc = 6

H  = 24
B  = 24

section_id = 1
from libdenavit.section.database import wide_flange_database
from libdenavit.section import SRC
d  = wide_flange_database[shape]['d']
bf = wide_flange_database[shape]['bf']
tf = wide_flange_database[shape]['tf']
tw = wide_flange_database[shape]['tw']

rho_sr = 0.01        
nbH = 2
nbB = 2
num_bars = 2*(nbH + nbB) - 4
Ab = rho_sr * H * B / num_bars
Dp = 3

section = SRC(H, B, d, bf, tf, tw, Fy, fc, 'US', nbH, nbB, Fylr, Dp, Ab=Ab)
steel_mat_type = 'ElasticPP'
reinf_mat_type = 'ElasticPP'
conc_mat_type = 'Concrete01_no_confinement'
nfy = 30
nfx = 30
#section.build_ops_fiber_section(section_id, steel_mat_type, reinf_mat_type, conc_mat_type, nfy, nfx)  
        

section.plot_section()

col = CrossSection2d(section, axis='x')

steel_mat_type = 'ElasticPP'
reinf_mat_type = 'ElasticPP'
conc_mat_type = 'Concrete01_no_confinement'
nfy = 30
nfx = 30
section_args = (steel_mat_type, reinf_mat_type, conc_mat_type, nfy, nfx) 
section_kwargs = {}

result_ops = col.run_ops_interaction(section_args=section_args, section_kwargs=section_kwargs)
plt.plot(result_ops["M1"], result_ops["P"], 'o-', label='OPS')
plt.xlabel('Bending Moment (kip-in.)')
plt.ylabel('Axial Force (kips)')
#plt.xlim(0, )
#plt.ylim(0, )
plt.legend()

print(section.p0)
print(section.Fy * section.As + section.Fylr * section.Asr + 1.0 * section.fc * section.Ac)

plt.show()
