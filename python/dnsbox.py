import numpy as np
from numpy import pi, sin, cos, exp
from numpy.fft import fftn, ifftn
import pickle

# Parameters:
global Nx, Ny, Nz, Nt, alphax, alphay, alphaz,\
       nu, Q,\
       iSaveRate1, iSaveRate2,\
       Deltak, kzero, uzero, CourantMin, CourantMax, t, tStepMax, tStepFix,\
       analytic, spherical, random

# Simulation parameters:
Nx = 64            # Spatial discretization
Ny = 64            # Spatial discretization
Nz = 64            # Spatial discretization
iSaveRate1 = 10    # Save state files every iSaveRate1 time steps
iSaveRate2 = 1     # Save state files every iSaveRate2 time steps
alphax = 1.0       # Base wave numbers
alphay = 1.0       # Base wave numbers
alphaz = 1.0       # Base wave numbers
nu = 0.2           # Viscosity           0.2
Q = 0.5            # Forcing amplitude   0.5
Deltak = 1.0       # k-window for shell averaging
CourantMin = 0.15  # Minimum Courant number
CourantMax = 0.2   # Maximum Courant number
tStepMax = 0.01    # Maximum time step

# Logicals
spherical = True   # Spherical truncation
analytic = False    # Analytical initial condition
random = True      # Random initial condition
tStepFix = True    # Fixed time-step simulation, if true

# Random initial field generation parameters:
kzero = 5.0
uzero = 1.0e-3


def setGrid():
    """
    Set grid in configuration and Fourier spaces
    """
    global Lx, Ly, Lz, x, y, z, kx, ky, kz
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
           intFact

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
    temphat = np.zeros((Nx, Ny, Nz), dtype='complex')
    uhattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    vhattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    whattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    utemp = np.zeros((Nx, Ny, Nz), dtype='float')
    vtemp = np.zeros((Nx, Ny, Nz), dtype='float')
    wtemp = np.zeros((Nx, Ny, Nz), dtype='float')

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
    global t, lamb, sk, sl, sm, A, kxm, kym, kzm, xx, yy, zz,\
           u, v, w, uhat, vhat, what, saveCount
    print('Initiating variables, grids, and fields')
    t = 0.0
    setGrid()
    alloc()
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                kxm[i, j, k] = kx[i]                                  # lint:ok
                kym[i, j, k] = ky[j]                                  # lint:ok
                kzm[i, j, k] = kz[k]                                  # lint:ok

                xx[i, j, k] = x[i]                                    # lint:ok
                yy[i, j, k] = y[j]                                    # lint:ok
                zz[i, j, k] = z[k]                                    # lint:ok
                if kx[i] == 0 and ky[j] == 0 and kz[k] == 0:
                    kkm[i, j, k] = 1.0e-13  # to avoid div by zeros
                else:
                    kkm[i, j, k] = np.real(kx[i] ** 2) \
                                 + np.real(ky[j] ** 2) \
                                 + np.real(kz[k] ** 2)  # k^2

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

        uhat[0, 0, 0] = 0.0
        vhat[0, 0, 0] = 0.0
        what[0, 0, 0] = 0.0

        project()
        uhat2u()

    else:
        print('read initial field from file')
        print('not implemented yet')
        return

    u2uhat()  # Fourier transform initial condition
    saveCount = 0

    return


def setTimeStep(tStep=tStepMax):
    global dt, intFact
    if tStepFix:
        dt = tStep
    else:
        print('Resetting time step')
        print('not implemented yet')
        return
    print('Set the time step to, dt = ', dt)
    #intFact = exp((- nu * (k2xm + k2ym + k2zm) + Q) * dt)
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
    global nonlinuhat, nonlinvhat, nonlinwhat, phat

    derivatives()  # compute derivatives
    utemp = np.real(ifftn(uhattemp))
    vtemp = np.real(ifftn(vhattemp))
    wtemp = np.real(ifftn(whattemp))

    # Compute N-S nonlinear terms:
    temp = utemp * ux + vtemp * uy + wtemp * uz
    nonlinuhat = fftn(temp)

    temp = utemp * vx + vtemp * vy + wtemp * vz
    nonlinvhat = fftn(temp)

    temp = utemp * wx + vtemp * wy + wtemp * wz
    nonlinwhat = fftn(temp)

    # Pressure
    phat = 1j * (kxm * nonlinuhat + kym * nonlinvhat + kzm * nonlinwhat) / kkm

    nonlinuhat = - nonlinuhat - 1j * kxm * phat
    nonlinvhat = - nonlinvhat - 1j * kym * phat
    nonlinwhat = - nonlinwhat - 1j * kzm * phat

    return


def timeStep(Nt):
    global u, v, w, uhat, vhat, what, uhattemp, vhattemp, whattemp, t
    print('Time stepping')

    for n in range(Nt):

        uhattemp = np.copy(uhat)
        vhattemp = np.copy(vhat)
        whattemp = np.copy(what)

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

        t = t + dt

    uhat2u()
    return


def Ekinetic():
    """
    Compute kinetic energy
    """
    Ekin = np.real(np.sum(np.conjugate(uhat) * uhat
                        + np.conjugate(vhat) * vhat
                        + np.conjugate(what) * what)) / (Nx * Ny * Nz) ** 2
    return Ekin


def checkError():

    print('t = ', t)

    uErrorMax = np.max((- A / (sk ** 2 + sl ** 2)) \
                     * (lamb * sl * cos(sk * xx) * sin(sl * yy) * sin(sm * zz)
                       + sm * sk * sin(sk * xx) * cos(sl * yy) * cos(sm * zz)) \
                     * exp(-1.0 * (lamb ** 2) * nu * t) - u)          # lint:ok

    vErrorMax = np.max((A / (sk ** 2 + sl ** 2)) \
                     * (lamb * sk * sin(sk * xx) * cos(sl * yy) * sin(sm * zz)
                       - sm * sl * cos(sk * xx) * sin(sl * yy) * cos(sm * zz)) \
                     * exp(-1.0 * (lamb ** 2) * nu * t) - v)          # lint:ok

    wErrorMax = np.max(A * cos(sk * xx) * cos(sl * yy) * sin(sm * zz) \
                     * exp(-1.0 * (lamb ** 2) * nu * t) - w)          # lint:ok

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


def saveStats():
    """
    Save flow statistics
    """
    k = Ekinetic()
    D = Dissipation()
    P = 2 * Q * Ekin
    C = Courant()
    div = divMax()

    f=open('stats.dat' ,'ab')
    f.write('%11.9f %11.9f %11.9f %11.9f %11.9f %11.9f\n' %
            (db.t,  k,      D,     P,     C,     div))
    f.close()

    return