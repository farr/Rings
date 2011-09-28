import numpy as np
import matplotlib.pyplot as pp

pp.rc('text', usetex=True)

data = np.loadtxt('/tmp/nature.dat')

ei = 3
ii = 4

ts = data[0::2, 0]
e1s = data[0::2, ei]
e2s = data[1::2, ei]
i1s = data[0::2, ii]
i2s = data[1::2, ii]

L1s = np.sqrt(1 - e1s**2.0)*np.cos(np.deg2rad(i1s))
L2s = np.sqrt(1 - e2s**2.0)*np.cos(np.deg2rad(i2s))

pp.subplot(4, 1, 1)
pp.subplots_adjust(hspace=0)
pp.plot(ts, i1s)
pp.axhline(90.0)
pp.gca().set_xticklabels([], visible=False)
pp.ylabel(r'$i_1$')

pp.subplot(4,1,2)
pp.yscale('log')
pp.plot(ts, 1 - e1s)
pp.gca().set_xticklabels([], visible=False)
pp.ylabel(r'$1-e_1$')

pp.subplot(4,1,3)
pp.plot(ts, L1s)
pp.gca().set_xticklabels([], visible=False)
pp.ylabel(r'$L_1^z$')

pp.subplot(4,1,4)
pp.plot(ts, L2s)
pp.ylabel(r'$L_2^z$')
