import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

runDir = "../test/0019"

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

ekin = np.loadtxt(runDir+"/Ekin.dat")

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
plt.figure(figsize=(8, 8))
plt.xlabel('$\mathcal{P}$', fontsize=24)
plt.ylabel('$\epsilon$', fontsize=24)
plt.plot(Production, Dissipation)

# Taylor microscale Reynolds number
plt.figure(figsize=(8,8))
plt.xlabel('$t$', fontsize=24)
plt.ylabel('$Re_{\lambda}$', fontsize=24)
plt.plot(time, Relambda)

# Turbulent reynolds number
plt.figure(figsize=(8,8))
plt.xlabel('$t$', fontsize=24)
plt.ylabel('$Re$', fontsize=24)
plt.plot(time, Returb)

plt.show()
