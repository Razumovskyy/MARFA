import numpy as np
import matplotlib.pyplot as plt

wn1, ka1 = np.loadtxt(fname='./calc/simpleSpectra100-110.dat', unpack=True, delimiter=',')
wn2, ka2 = np.loadtxt(fname='./calc/simpleSpectraCutOff100-110.dat', unpack=True, delimiter=',')


fig, ax = plt.subplots()

ax.plot(wn1, ka1, color='r', label='only inbound spectral lines')
ax.plot(wn2, ka2, color='b', label='outbound with cut-off 25 cm-1')

ax.set_xlabel(r'wavenumber, $\mathregular{cm^{-1}}$')
ax.set_ylabel(r'$k,\,\mathregular{km^{-1}}$')
ax.grid(which='major', axis='both', color='gray', alpha=0.5)

plt.title(r'H$_2$O Absorption Coefficient in $\mathregular{km^{-1}}$ at 1 atm in the 2100-2110 $\mathregular{cm^{-1}}$ Range for simple Lorentz shape')
plt.tight_layout()
plt.show()