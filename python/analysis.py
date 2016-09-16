import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

runDir = "../test/0106"
Ni = 0
Nf = -1

#Read simulation parameters:
with open(runDir+"/info.dat") as fp:
    for i, line in enumerate(fp):
        if i == 3:
            alphax = float(line.split()[2])
            alphay = float(line.split()[3])
            alphaz = float(line.split()[4])
        if i == 5:
            nu = float(line.split()[2])
        if i == 6:
            Q = float(line.split()[2])

ekin = np.loadtxt(runDir+"/Ekin.dat")[Ni:Nf]

#Quantities:
time = ekin[:, 0]
KineticEnergy = ekin[:, 1]
Dissipation = ekin[:, 2]
Production = ekin[:, 3]
urms = np.sqrt(2.0 * KineticEnergy / 3.0)
lint = urms ** 3 / Dissipation # Integral length scale:
lambdaTaylor = np.sqrt(15.0 * nu * (Dissipation ** (-1.0))) * urms 
Relambda = lambdaTaylor * urms / nu 
Returb = lint * urms / nu 

 
# Turbulent kinetic energy
plt.figure(figsize=(8, 8))
plt.xlabel('$t$', fontsize=24)
plt.ylabel('$k$', fontsize=24)
plt.plot(time, KineticEnergy)

# Production-Dissipation plot
plt.figure(figsize=(6, 6))
plt.xlabel('$\mathcal{P}$', fontsize=36)
plt.ylabel('$\epsilon$', fontsize=36)
plt.plot(Production, Dissipation)
ax = plt.gca()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(24)
    
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(24)

# Taylor microscale Reynolds number
plt.figure(figsize=(8,8))
plt.xlabel('$t$', fontsize=16)
plt.ylabel('$Re_{\lambda}$', fontsize=16)
plt.plot(time, Relambda)

# Turbulent reynolds number
plt.figure(figsize=(8,8))
plt.xlabel('$t$', fontsize=24)
plt.ylabel('$Re$', fontsize=24)
plt.plot(time, Returb)

plt.show()
