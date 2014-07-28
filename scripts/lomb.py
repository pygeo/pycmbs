from pycmbs.statistic import lomb_scargle_periodogram
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')


t = np.arange(1, 365*10, 1.)
y = np.random.random(len(t)) * 5.
y = 5. * np.cos(2.*np.pi * t / 365)

P = np.linspace(0.001, 1000., 500)

A, B = lomb_scargle_periodogram(t, P, y)

f = plt.figure()
ax1 = f.add_subplot(211)
ax2 = f.add_subplot(212)

ax1.plot(P, A)
ax2.plot(P, B)
ax2.set_xlabel('period [days]')
ax2.set_ylabel('phase [rad]')
ax1.set_ylabel('amplitude')
ax1.grid()
ax2.grid()

plt.show()
