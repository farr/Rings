import numpy as np
import matplotlib.pyplot as pp
import scipy as scp

pp.rc('text', usetex=True)

data = np.loadtxt('/tmp/nature.dat')

mi = 1
ai = 2
ei = 3
ii = 4
Oi = 5
oi = 6
si = 11
Ii = 9
Ri = 10

def lorb(data):
    m = data[mi]
    a = data[ai]
    e = data[ei]
    i = np.pi/180.0*data[ii]
    O = np.pi/180.0*data[Oi]

    mm=np.sqrt((1.0+m)/(a*a*a))

    Lmag=m*mm*a*a*np.sqrt(1.0 - e*e)

    return Lmag*np.array([np.cos(O-np.pi/2)*np.sin(i), 
                          np.sin(O-np.pi/2)*np.sin(i),
                          np.cos(i)])

def lspin(data):
    s = data[si:(si+3)]
    I = data[Ii]

    return I*s

def ltot(data):
    return lspin(data[0,:]) + lspin(data[1,:]) + lspin(data[2,:]) + lorb(data[1,:]) + lorb(data[2,:])

nrow=4
ncol=1
iplot=1

ts = data[::3, 0]
e1s = data[1::3, ei]
e2s = data[2::3, ei]
i1s = data[1::3, ii]
i2s = data[2::3, ii]
a1s = data[1::3, ai]
a2s = data[2::3, ai]
s0s = data[0::3, si:(si+3)]
s1s = data[1::3, si:(si+3)]

thetas=[]
for i in range(s0s.shape[0]):
    s0=s0s[i, :]
    s1=s1s[i, :]

    s0 /= np.linalg.norm(s0)
    s1 /= np.linalg.norm(s1)

    thetas.append(180.0/np.pi*np.arccos(np.dot(s0,s1)))
thetas=np.array(thetas)

tis=np.linspace(min(ts), max(ts), 1000)

pp.subplot(nrow,ncol,iplot)
iplot += 1
pp.subplots_adjust(hspace=0)
pp.plot(tis, scp.interp(tis, ts, i1s), 'k-', tis, scp.interp(tis, ts, thetas), 'r-')
pp.axhline(90.0)
pp.gca().set_xticklabels([], visible=False)
pp.ylabel(r'$i_1$')

pp.subplot(nrow,ncol,iplot)
iplot += 1
pp.yscale('log')
pp.plot(tis, scp.interp(tis, ts, 1 - e1s))
pp.gca().set_xticklabels([], visible=False)
pp.ylabel(r'$1-e_1$')

pp.subplot(nrow,ncol,iplot)
iplot += 1
pp.plot(tis, scp.interp(tis, ts, a1s*(1-e1s)), 'r-', tis, scp.interp(tis, ts, a2s*(1-e2s)), 'g-')
pp.yscale('log')
pp.ylabel(r'$r_p$')

pp.subplot(nrow,ncol,iplot)
l0=ltot(data[0:3, :])
l0norm=np.linalg.norm(l0)
ls=np.array([np.linalg.norm(ltot(data[3*i:(3*i+3),:])-l0)/l0norm for i in range(data.shape[0]/3)])
pp.plot(tis, scp.interp(tis, ts, ls), 'k')
pp.yscale('log')
pp.ylabel(r'$\left| \vec{L}_0^\mathrm{tot} - \vec{L}^\mathrm{tot} \right|/\left|\vec{L}_0^\mathrm{tot}\right|$')
pp.xlabel(r'$t$')

pp.savefig('/tmp/nature.pdf')

pp.show()
