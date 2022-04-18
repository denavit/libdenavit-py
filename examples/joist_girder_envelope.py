import numpy as np
import matplotlib.pyplot as plt
from libdenavit import JoistGirder

d_in = 36
L_ft = 40
num_spaces = 8
Pa_kips = 20
joist_girder = JoistGirder('LRFD',L_ft,d_in,num_spaces,Pa_kips)

x = np.zeros(10*num_spaces)
for i in range(num_spaces):
    x[(10*i):(10*(i+1))] = L_ft*np.linspace(i/num_spaces, (i+1)/num_spaces, num=10)

(pos_env_M,neg_env_M) = joist_girder.moment_strength_envelope(x,units='kip-ft')
(pos_env_V,neg_env_V) = joist_girder.shear_strength_envelope(x,units='kips')

plt.figure()
plt.plot(x,pos_env_M)
plt.plot(x,neg_env_M)
plt.xlabel('Position along length of joist girder (ft)')
plt.ylabel('Moment (kip-ft)')
plt.xlim(0,L_ft)

plt.figure()
plt.plot(x,pos_env_V)
plt.plot(x,neg_env_V)
plt.xlabel('Position along length of joist girder (ft)')
plt.ylabel('Shear (kips)')
plt.xlim(0,L_ft)

plt.show()
