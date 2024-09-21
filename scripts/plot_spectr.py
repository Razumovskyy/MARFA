import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'

wn1, ka1 = np.loadtxt(fname='SPECTR', skiprows=1, unpack=True)
wn2, ka2 = np.loadtxt(fname='SPECTR_PYTHON.dat', skiprows=1, unpack=True)


fig, ax = plt.subplots()

# ax.plot(wn1, ka1, color='b', linestyle='-', marker='o', markersize=10, label='Voigt (125)')
ax.plot(wn1, ka1, color='b', linestyle='-', label='Voigt (125)')

# ax.plot(wn2, ka2, color='g', linestyle='-', marker='s', markersize=10, label='Tonkov (350)')
ax.plot(wn2, ka2, color='g', linestyle='-', label='Tonkov (350)')

ax.set_xlabel(r'Wavenumber [$\mathregular{cm^{-1}}$]')
ax.set_ylabel(r'Absorption cross-section [cm$^2$ mol$^{-1}$]')
ax.grid(which='major', axis='both', color='gray', alpha=0.5)

# Configure major and minor ticks
ax.minorticks_on()  # Add minor ticks
ax.tick_params(axis='y', which='both', right=True)  # Tick marks on the right

# ax.legend()
# ax.set_xlim(left=7500, right=10000)

plt.tight_layout()
plt.show()
