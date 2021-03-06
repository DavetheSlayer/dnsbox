import numpy as np
from numpy import pi, sin, cos, exp
from numpy.fft import fftn, ifftn
import pickle
from pylab import *

# Parameters:
global Nx, Ny, Nz, Nt, alphax, alphay, alphaz,\
       nu, Q,\
       iSaveRate1, iSaveRate2,\
       Deltak, kzero, uzero, CourantMin, CourantMax, t, tStepMax, tStepFix,\
       analytic, spherical, random, dt

# Simulation parameters:
Nx = 64            # Spatial discretization
Ny = 64            # Spatial discretization
Nz = 64            # Spatial discretization
iSaveRate1 = 10    # Save state files every iSaveRate1 time steps
iSaveRate2 = 1     # Save state files every iSaveRate2 time steps
alphax = 1.0       # Base wave numbers
alphay = 1.0       # Base wave numbers
alphaz = 1.0       # Base wave numbers
nu = 0.07          # Viscosity           0.2
Q = 0.0            # Forcing amplitude   0.5
Deltak = 1.0       # k-window for shell averaging
CourantMin = 0.15  # Minimum Courant number
CourantMax = 0.2   # Maximum Courant number
tStepMax = 0.01    # Maximum time step
bandlim = True     # Band-limited forcing if true
kF = 2.5           # wave number cut off for band limited forcing
epsW = 0.1         # Rate of energy input for band limited forcing

# Logicals
spherical = False  # Spherical truncation
analytic = False   # Analytical initial condition
random = True      # Random initial condition
tStepFix = False   # Fixed time-step simulation, if true

# Random initial field generation parameters:
kzero = 3.0
uzero = 1.0e0


def setGrid():
    """
    Set grid in configuration and Fourier spaces
    """
    global Lx, Ly, Lz, x, y, z, kx, ky, kz, kspec, Espec
    Lx = 2 * pi * alphax
    Ly = 2 * pi * alphay
    Lz = 2 * pi * alphaz

    x = np.array([i * (Lx / Nx) for i in range(-Nx / 2, 1 + Nx / 2)],
                 dtype='float')
    y = np.array([i * (Ly / Ny) for i in range(-Ny / 2, 1 + Ny / 2)],
                 dtype='float')
    z = np.array([i * (Lz / Nz) for i in range(-Nz / 2, 1 + Nz / 2)],
                 dtype='float')
    kx = alphax * np.array([n for n in list(range(0, Nx / 2)) + [0]
                                     + list(range(-Nx / 2 + 1, 0))])
    ky = alphay * np.array([n for n in list(range(0, Ny / 2)) + [0]
                                     + list(range(-Ny / 2 + 1, 0))])
    kz = alphaz * np.array([n for n in list(range(0, Nz / 2)) + [0]
                                     + list(range(-Nz / 2 + 1, 0))])


    kspec = np.arange(Deltak, np.max(kx) * (2.0/3.0) + Deltak, Deltak)
    Espec = np.zeros(kspec.shape)


    return


def alloc():
    """
    Allocate arrays
    """
    global u, v, w, ux, uy, uz, vx, vy, vz, wx, wy, wz,\
           omegax, omegay, omegaz,\
           uhat, vhat, what, phat, nonlinuhat, nonlinvhat, nonlinwhat,\
           xx, yy, zz, kxm, kym, kzm, kkm, \
           temp, temphat, utemp, vtemp, wtemp, uhattemp, vhattemp, whattemp,\
           intFact, dealias, force

    u = np.zeros((Nx, Ny, Nz), dtype='float')
    v = np.zeros((Nx, Ny, Nz), dtype='float')
    w = np.zeros((Nx, Ny, Nz), dtype='float')

    ux = np.zeros((Nx, Ny, Nz), dtype='float')
    uy = np.zeros((Nx, Ny, Nz), dtype='float')
    uz = np.zeros((Nx, Ny, Nz), dtype='float')
    vx = np.zeros((Nx, Ny, Nz), dtype='float')
    vy = np.zeros((Nx, Ny, Nz), dtype='float')
    vz = np.zeros((Nx, Ny, Nz), dtype='float')
    wx = np.zeros((Nx, Ny, Nz), dtype='float')
    wy = np.zeros((Nx, Ny, Nz), dtype='float')
    wz = np.zeros((Nx, Ny, Nz), dtype='float')

    omegax = np.zeros((Nx, Ny, Nz), dtype='float')
    omegay = np.zeros((Nx, Ny, Nz), dtype='float')
    omegaz = np.zeros((Nx, Ny, Nz), dtype='float')

    uhat = np.zeros((Nx, Ny, Nz), dtype='complex')
    vhat = np.zeros((Nx, Ny, Nz), dtype='complex')
    what = np.zeros((Nx, Ny, Nz), dtype='complex')
    phat = np.zeros((Nx, Ny, Nz), dtype='complex')

    nonlinuhat = np.zeros((Nx, Ny, Nz), dtype='complex')
    nonlinvhat = np.zeros((Nx, Ny, Nz), dtype='complex')
    nonlinwhat = np.zeros((Nx, Ny, Nz), dtype='complex')

    xx = np.zeros((Nx, Ny, Nz), dtype='float')
    yy = np.zeros((Nx, Ny, Nz), dtype='float')
    zz = np.zeros((Nx, Ny, Nz), dtype='float')

    kxm = np.zeros((Nx, Ny, Nz), dtype='complex')
    kym = np.zeros((Nx, Ny, Nz), dtype='complex')
    kzm = np.zeros((Nx, Ny, Nz), dtype='complex')

    kkm = np.zeros((Nx, Ny, Nz), dtype='float')

    intFact = np.zeros((Nx, Ny, Nz), dtype='float')
    temp = np.zeros((Nx, Ny, Nz), dtype='float')
    temphat = np.zeros((Nx, Ny, Nz), dtype='complex')
    uhattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    vhattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    whattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    utemp = np.zeros((Nx, Ny, Nz), dtype='float')
    vtemp = np.zeros((Nx, Ny, Nz), dtype='float')
    wtemp = np.zeros((Nx, Ny, Nz), dtype='float')

    dealias = np.ones((Nx, Ny, Nz), dtype='float')
    force = np.ones((Nx, Ny, Nz), dtype='float')

    return


def u2uhat():
    """
    Forward Fourier transform the fields
    """
    global u, v, w, uhat, vhat, what
    uhat = fftn(u)                                                    # lint:ok
    vhat = fftn(v)                                                    # lint:ok
    what = fftn(w)                                                    # lint:ok
    return


def uhat2u():
    """
    Forward Fourier transform the fields
    """
    global u, v, w, uhat, vhat, what
    u = np.real(ifftn(uhat))                                          # lint:ok
    v = np.real(ifftn(vhat))                                          # lint:ok
    w = np.real(ifftn(what))                                          # lint:ok
    return


def project():
    """
    Project state to make it divergence free
    """
    global uhat, vhat, what, uhattemp, vhattemp, whattemp
    uhattemp = np.copy(uhat)
    vhattemp = np.copy(vhat)
    whattemp = np.copy(what)

    uhat = uhattemp - kxm * kxm * uhattemp / kkm \
                    - kxm * kym * vhattemp / kkm \
                    - kxm * kzm * whattemp / kkm

    vhat = vhattemp - kym * kxm * uhattemp / kkm \
                    - kym * kym * vhattemp / kkm \
                    - kym * kzm * whattemp / kkm

    what = whattemp - kzm * kxm * uhattemp / kkm \
                    - kzm * kym * vhattemp / kkm \
                    - kzm * kzm * whattemp / kkm

    return


def init():
    """
    Initiate simulation
    """
    global t, lamb, sk, sl, sm, A, kxm, kym, kzm, xx, yy, zz, \
           u, v, w, uhat, vhat, what, saveCount, dt, intFact, \
           uhattemp, vhattemp, whattemp, dealias, force
    print('Initiating variables, grids, and fields')
    t = 0.0
    setGrid()
    alloc()
    kxMax = np.max(kx)
    kyMax = np.max(ky)
    kzMax = np.max(kz)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                kxm[i, j, k] = kx[i]                                  # lint:ok
                kym[i, j, k] = ky[j]                                  # lint:ok
                kzm[i, j, k] = kz[k]                                  # lint:ok

                xx[i, j, k] = x[i]                                    # lint:ok
                yy[i, j, k] = y[j]                                    # lint:ok
                zz[i, j, k] = z[k]                                    # lint:ok

                kk = np.real(kx[i] ** 2) \
                   + np.real(ky[j] ** 2) \
                   + np.real(kz[k] ** 2)  # k^2

                if kx[i] == 0 and ky[j] == 0 and kz[k] == 0:
                    kkm[i, j, k] = 1.0e-13  # to avoid div by zeros
                else:
                    kkm[i, j, k] = kk

                # dealias array sets |kk| = 0 and |kk| > |kkmax| modes to zero
                # when element-wise multiplied an Fourier-space array.
                if spherical and (np.sqrt(kk) >= (2.0 / 3.0) * kxMax
                               or kk == 0.0):

                    dealias[i, j, k] = 0.0

                elif (kx[i] >= (2.0 / 3.0) * kxMax or
                      ky[j] >= (2.0 / 3.0) * kyMax or
                      kz[k] >= (2.0 / 3.0) * kzMax or
                      kk == 0.0):

                    dealias[i, j, k] = 0.0

                if np.sqrt(kk) >= kF:

                    force[i, j, k] = 0.0


    if analytic:
        # Initiate simulation from an analytical solution for test. See
        # https://en.wikibooks.org/wiki/Parallel_Spectral_Numerical_Methods/
        # The_Two-_and_Three-Dimensional_Navier-Stokes_Equations
        sk = 1.0
        sl = 1.0
        sm = 1.0
        lamb = np.sqrt(sk ** 2 + sl ** 2 + sm ** 2)
        A = 1.0
        print('generating initial condition based on analytical solution')
        u = (- A / (sk ** 2 + sl ** 2)) \
          * (lamb * sl * cos(sk * xx) * sin(sl * yy) * sin(sm * zz)
           + sm * sk * sin(sk * xx) * cos(sl * yy) * cos(sm * zz)) \
          * exp(-1.0 * (lamb ** 2) * nu * t)                          # lint:ok
        v = (A / (sk ** 2 + sl ** 2)) \
          * (lamb * sk * sin(sk * xx) * cos(sl * yy) * sin(sm * zz)
           - sm * sl * cos(sk * xx) * sin(sl * yy) * cos(sm * zz)) \
          * exp(-1.0 * (lamb ** 2) * nu * t)                          # lint:ok
        w = A * cos(sk * xx) * cos(sl * yy) * sin(sm * zz) \
          * exp(-1.0 * (lamb ** 2) * nu * t)                          # lint:ok
        # del xx, yy, zz                                              # lint:ok
    elif random:
        # Initiation method: Rosales & Meneveau (2005), Phys. Fluids. eq.9
        print('random initial field')
        uhat = np.sqrt(4 * (uzero ** 2) * kkm * (kzero ** (-5))
                     * (pi ** (-3.0 / 2.0)) * (2.0 ** (-0.5))
                     * exp(-2.0 * kkm / (kzero ** 2))) \
               * exp(1j * (np.random.rand(Nx, Ny, Nz) - 0.5) * 2 * pi) \
               * Nx * Ny * Nz
        vhat = np.sqrt(4 * (uzero ** 2) * kkm * (kzero ** (-5))
                     * (pi ** (-3.0 / 2.0)) * (2.0 ** (-0.5))
                     * exp(-2.0 * kkm / (kzero ** 2))) \
               * exp(1j * (np.random.rand(Nx, Ny, Nz) - 0.5) * 2 * pi) \
               * Nx * Ny * Nz
        what = np.sqrt(4 * (uzero ** 2) * kkm * (kzero ** (-5))
                     * (pi ** (-3.0 / 2.0)) * (2.0 ** (-0.5))
                     * exp(-2.0 * kkm / (kzero ** 2))) \
               * exp(1j * (np.random.rand(Nx, Ny, Nz) - 0.5) * 2 * pi) \
               * Nx * Ny * Nz

        uhat = dealias * uhat
        vhat = dealias * vhat
        what = dealias * what
        project()
        uhat2u()

    else:
        print('read initial field from file')
        loadState('state0000.pkl')

    dt = tStepMax  # Initial time-step setting
    intFact = exp((- nu * kkm + Q) * dt)  # Set the integrating factor
    u2uhat()  # Fourier transform initial condition
    uhat = dealias * uhat
    vhat = dealias * vhat
    what = dealias * what
    uhattemp = np.copy(uhat)
    vhattemp = np.copy(vhat)
    whattemp = np.copy(what)
    saveCount = 0

    return


def setTimeStep(tStep=tStepMax):
    """
    Change the time-step such that the Courant number is set to
    (CourantMin+CourantMax)/2
    This function should be called after the Courant number is computed when
    running in adaptive time-stepping mode
    """
    global dt, intFact
    if tStepFix:
        dt = tStep
    elif CourantNumber > CourantMax or (CourantNumber < CourantMin
                                        and dt < tStepMax):
        dt = ((CourantMax + CourantMin) / (2.0 * CourantNumber)) * dt
        if (dt > tStepMax):
            dt = tStepMax
    else:
        return

    print('Set the time step to, dt = ', dt)
    intFact = exp((- nu * kkm + Q) * dt)
    return


def derivatives():
    """
    Compute space derivatives from (u,v,w)hattemp
    """
    global ux, uy, uz, vx, vy, vz, wx, wy, wz

    temphat = 1j * kxm * uhattemp
    ux = np.real(ifftn(temphat))
    temphat = 1j * kym * uhattemp
    uy = np.real(ifftn(temphat))
    temphat = 1j * kzm * uhattemp
    uz = np.real(ifftn(temphat))

    temphat = 1j * kxm * vhattemp
    vx = np.real(ifftn(temphat))
    temphat = 1j * kym * vhattemp
    vy = np.real(ifftn(temphat))
    temphat = 1j * kzm * vhattemp
    vz = np.real(ifftn(temphat))

    temphat = 1j * kxm * whattemp
    wx = np.real(ifftn(temphat))
    temphat = 1j * kym * whattemp
    wy = np.real(ifftn(temphat))
    temphat = 1j * kzm * whattemp
    wz = np.real(ifftn(temphat))

    return


def nonLinear():
    """
    Compute nonlinear term for (u,v,w)hattemp
    """
    global nonlinuhat, nonlinvhat, nonlinwhat, phat, omegax, omegay, omegaz, \
           utemp, vtemp, wtemp, temp

    derivatives()  # compute derivatives
    utemp = np.real(ifftn(uhattemp))
    vtemp = np.real(ifftn(vhattemp))
    wtemp = np.real(ifftn(whattemp))
    
    omegax = wy - vz
    omegay = uz - wx
    omegaz = vx - uy
    
    # Compute N-S nonlinear terms:
    # temp = utemp * ux + vtemp * uy + wtemp * uz
    temp = omegay * wtemp - omegaz * vtemp
    nonlinuhat = fftn(temp)
    
    # temp = utemp * vx + vtemp * vy + wtemp * vz
    temp = omegaz * utemp - omegax * wtemp
    nonlinvhat = fftn(temp)
    
    # temp = utemp * wx + vtemp * wy + wtemp * wz
    temp = omegax * vtemp - omegay * utemp
    nonlinwhat = fftn(temp)

    # Pressure
    phat = 1j * (kxm * nonlinuhat + kym * nonlinvhat + kzm * nonlinwhat) / kkm

    if bandlim:
        Eband = 0.5 * np.real(np.sum(np.conjugate(uhattemp) * (force * uhattemp)
                                   + np.conjugate(vhattemp) * (force * vhattemp)
                                   + np.conjugate(whattemp) * (force * whattemp))
                                    / (Nx * Ny * Nz) ** 2)

        nonlinuhat = dealias * (- nonlinuhat - 1j * kxm * phat 
                               + (epsW / (2.0 * Eband)) * uhattemp)
        nonlinvhat = dealias * (- nonlinvhat - 1j * kym * phat
                               + (epsW / (2.0 * Eband)) * vhattemp)
        nonlinwhat = dealias * (- nonlinwhat - 1j * kzm * phat
                               + (epsW / (2.0 * Eband)) * whattemp)

    else:
        nonlinuhat = dealias * (- nonlinuhat - 1j * kxm * phat)
        nonlinvhat = dealias * (- nonlinvhat - 1j * kym * phat)
        nonlinwhat = dealias * (- nonlinwhat - 1j * kzm * phat)

    return


def timeStep(Nt):
    global u, v, w, uhat, vhat, what, uhattemp, vhattemp, whattemp, t
    print('Time stepping')

    for n in range(Nt):

        nonLinear()  # Compute nonlinear term

        #Predictor-Corrector:

        uhattemp = intFact * (uhat + dt * nonlinuhat)   # Prediction
        uhat = intFact * (uhat + dt * nonlinuhat * 0.5) # contribution from
                                                        # prediction
        vhattemp = intFact * (vhat + dt * nonlinvhat)   # Prediction
        vhat = intFact * (vhat + dt * nonlinvhat * 0.5) # contribution from
                                                        # prediction
        whattemp = intFact * (what + dt * nonlinwhat)   # Prediction
        what = intFact * (what + dt * nonlinwhat * 0.5) # contribution from
                                                        # prediction

        nonLinear()  # Compute nonlinear term again

        uhat = uhat + dt * nonlinuhat * 0.5
        vhat = vhat + dt * nonlinvhat * 0.5
        what = what + dt * nonlinwhat * 0.5

        # project()  # Apply projection operator to ensure div-free

        uhattemp = np.copy(uhat)
        vhattemp = np.copy(vhat)
        whattemp = np.copy(what)

        t = t + dt

    uhat2u()
    return


def checkError():

    print('t = ', t)

    uErrorMax = np.max((- A / (sk ** 2 + sl ** 2))
                     * (lamb * sl * cos(sk * xx) * sin(sl * yy) * sin(sm * zz)
                       + sm * sk * sin(sk * xx) * cos(sl * yy) * cos(sm * zz))
                     * exp(-1.0 * (lamb ** 2) * nu * t) - u)

    vErrorMax = np.max((A / (sk ** 2 + sl ** 2))
                     * (lamb * sk * sin(sk * xx) * cos(sl * yy) * sin(sm * zz)
                       - sm * sl * cos(sk * xx) * sin(sl * yy) * cos(sm * zz))
                     * exp(-1.0 * (lamb ** 2) * nu * t) - v)

    wErrorMax = np.max(A * cos(sk * xx) * cos(sl * yy) * sin(sm * zz)
                     * exp(-1.0 * (lamb ** 2) * nu * t) - w)

    print('Maximum errors at the final time step are: ')
    print(uErrorMax, vErrorMax, wErrorMax)
    return


def saveState():

    global saveCount
    f = open('state'+str(saveCount).zfill(4)+'.pkl', 'wb')
    pickle.dump([u, v, w], f)
    f.close()
    saveCount = saveCount + 1
    return


def loadState(stateName):

    global u, v, w
    f = open(stateName, 'rb')
    u, v, w = pickle.load(f)
    f.close()
    return


def stats():
    """
    Save flow statistics
    """
    global CourantNumber
    derivatives()

    CourantNumber = np.max(np.abs(u)) * dt / (Lx / Nx) \
                  + np.max(np.abs(v)) * dt / (Ly / Ny) \
                  + np.max(np.abs(w)) * dt / (Lz / Nz)

    div = np.max(ux + vy + wz)  # Maximum divergence

    k = 0.5 * np.real(np.sum(np.conjugate(uhat) * uhat
                           + np.conjugate(vhat) * vhat
                           + np.conjugate(what) * what)) / (Nx * Ny * Nz) ** 2

    D = np.sum(ux ** 2 + uy ** 2 + uz ** 2
             + vx ** 2 + vy ** 2 + vz ** 2
             + wx ** 2 + wy ** 2 + wz ** 2) * nu / (Nx * Ny * Nz)

    P = 2 * Q * k

    C = CourantNumber

    f = open('stats.dat', 'ab')
    f.write('%.18e %.18e %.18e %.18e %.18e %.18e\n' %
            (t,     k,     D,     P,     C,     div))
    #np.savetxt(f, np.array([t, k, D, P, C, div]), delimiter=' ')
    f.close()

    return

def spectrum(plot=True):

    global Espec, Ezero
    Ezero = 0.0
    Espec[:] = 0.0
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):

                absk = np.sqrt(np.real(kx[i] ** 2) \
                     + np.real(ky[j] ** 2) \
                     + np.real(kz[k] ** 2))  # |k|
                if absk < int(absk / Deltak) + 0.5 * Deltak:
                    nk = int(absk / Deltak) - 1.0  # python counts from 0
                else:
                    nk = int(absk / Deltak)

                if absk == 0.0:
                    Ezero = Ezero + 0.5 * np.real(np.sum(
                np.conjugate(uhat[i,j,k]) * uhat[i,j,k]
              + np.conjugate(vhat[i,j,k]) * vhat[i,j,k]
              + np.conjugate(what[i,j,k]) * what[i,j,k])) \
              / (Nx * Ny * Nz) ** 2

                if nk >= 0 and nk < len(kspec):
                    Espec[nk] = Espec[nk] + 0.5 * np.real(np.sum(
                np.conjugate(uhat[i,j,k]) * uhat[i,j,k]
              + np.conjugate(vhat[i,j,k]) * vhat[i,j,k]
              + np.conjugate(what[i,j,k]) * what[i,j,k])) \
              / (Nx * Ny * Nz) ** 2 * Deltak

    print("Ezero = ", Ezero)
    loglog(kspec, Espec)
    show()

    return

