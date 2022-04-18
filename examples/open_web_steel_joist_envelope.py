import numpy as np
import matplotlib.pyplot as plt
from libdenavit import OpenWebSteelJoist

L_ft  = 40
wTL_plf = 100
wLL_plf = 80
joist = OpenWebSteelJoist('LRFD',L_ft,wTL_plf,wLL_plf)
joist.minimum_shear_reversal_strength_ratio = 0.125

x = np.concatenate((np.linspace(0.0, 0.5*L_ft, num=21),np.linspace(0.5*L_ft, L_ft, num=21)))

(pos_env_M,neg_env_M) = joist.moment_strength_envelope(x)
(pos_env_V,neg_env_V) = joist.shear_strength_envelope(x)

plt.figure()
plt.plot(x,pos_env_M)
plt.plot(x,neg_env_M)
plt.xlabel('Position along length of joist (ft)')
plt.ylabel('Moment (lb-ft)')
plt.xlim(0,L_ft)

plt.figure()
plt.plot(x,pos_env_V)
plt.plot(x,neg_env_V)
plt.xlabel('Position along length of joist (ft)')
plt.ylabel('Shear (lbs)')
plt.xlim(0,L_ft)

plt.show()
