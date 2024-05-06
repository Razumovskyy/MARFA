import numpy as np
import matplotlib.pyplot as plt

wn, ka = np.loadtxt(fname='SPECTR', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(wn, ka, color='r')
ax.set_xlabel(r'wavenumber, $\mathregular{cm^{-1}}$')
ax.set_ylabel(r'$\log(\alpha),\,\mathregular{km^{-1}}$')
ax.grid(which='major', axis='both', color='gray', alpha=0.5)

plt.title(r'SO$_2$ Absorption Coefficient at 98 km in the 400-600 $\mathregular{cm^{-1}}$ Range (HITRAN, Voigt Shape)')
plt.tight_layout()
plt.show()