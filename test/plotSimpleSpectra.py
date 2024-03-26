import numpy as np
import matplotlib.pyplot as plt

wn, ka = np.loadtxt(fname='./simpleSpectra.dat', unpack=True, delimiter=',')

fig, ax = plt.subplots()

ax.plot(wn, ka, color='r')
ax.set_xlabel(r'wavenumber, $\mathregular{cm^{-1}}$')
ax.set_ylabel(r'$k,\,\mathregular{km^{-1}}$')
ax.grid(which='major', axis='both', color='gray', alpha=0.5)

plt.title(r'H$_2$O Absorption Coefficient in $\mathregular{km^{-1}}$ at 1 atm in the 100-110 $\mathregular{cm^{-1}}$ Range for simple Lorentz shape')
plt.tight_layout()
plt.show()