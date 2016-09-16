import numpy as np
from numpy import pi, sin, cos, exp
from numpy.fft import fftn, ifftn

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
Nt = 10            # Maximum number of time stepes
iSaveRate1 = 10    # Save state files every iSaveRate1 time steps
iSaveRate2 = 1     # Save state files every iSaveRate2 time steps
alphax = 1.0       # Base wave numbers
alphay = 1.0       # Base wave numbers
alphaz = 1.0       # Base wave numbers
nu = 1.0           # Viscosity           0.2
Q = 0.0            # Forcing amplitude   0.5
Deltak = 1.0       # k-window for shell averaging
CourantMin = 0.15  # Minimum Courant number
CourantMax = 0.2   # Maximum Courant number
tStepMax = 0.01    # Maximum time step

# Logicals
spherical = True   # Spherical truncation
analytic = True    # Analytical initial condition
random = False     # Random initial condition
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
    kx = alphax * np.array([1j * n for n in list(range(0, Nx / 2)) + [0]
                                          + list(range(-Nx / 2 + 1, 0))])
    ky = alphay * np.array([1j * n for n in list(range(0, Ny / 2)) + [0]
                                          + list(range(-Ny / 2 + 1, 0))])
    kz = alphaz * np.array([1j * n for n in list(range(0, Nz / 2)) + [0]
                                          + list(range(-Nz / 2 + 1, 0))])

    return


def alloc():
    """
    Allocate arrays
    """
    global u, v, w, ux, uy, uz, vx, vy, vz, wx, wy, wz,\
           omegax, omegay, omegaz,\
           uhat, vhat, what, phat, nonlinuhat, nonlinvhat, nonlinwhat,\
           xx, yy, zz, kxm, kym, kzm, k2xm, k2ym, k2zm,\
           temp, temphat, utemp, vtemp, wtemp, uhattemp, vhattemp, whattemp,\
           intFact  # , uhatold, vhatold, whatold

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

    k2xm = np.zeros((Nx, Ny, Nz), dtype='float')
    k2ym = np.zeros((Nx, Ny, Nz), dtype='float')
    k2zm = np.zeros((Nx, Ny, Nz), dtype='float')

    intFact = np.zeros((Nx, Ny, Nz), dtype='float')
    temphat = np.zeros((Nx, Ny, Nz), dtype='complex')
    uhattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    vhattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    whattemp = np.zeros((Nx, Ny, Nz), dtype='complex')
    #uhatold = np.zeros((Nx, Ny, Nz), dtype='complex')
    #vhatold = np.zeros((Nx, Ny, Nz), dtype='complex')
    #whatold = np.zeros((Nx, Ny, Nz), dtype='complex')
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


def init():
    """
    Initiate simulation
    """
    global t, lamb, sk, sl, sm, A, kxm, kym, kzm, k2xm, k2ym, k2zm, xx, yy, zz,\
           u, v, w, uhat, vhat, what
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

                k2xm[i, j, k] = np.real(kx[i] ** 2)                   # lint:ok
                k2ym[i, j, k] = np.real(ky[i] ** 2)                   # lint:ok
                k2zm[i, j, k] = np.real(kz[i] ** 2)                   # lint:ok

                xx[i, j, k] = x[i]                                    # lint:ok
                yy[i, j, k] = y[i]                                    # lint:ok
                zz[i, j, k] = z[i]                                    # lint:ok

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
        u = (-A / (sk ** 2 + sl ** 2)) \
          * (lamb * sl * cos(sk * xx) * sin(sl * yy) * sin(sm * zz)
           + sm * sk * sin(sk * xx) * cos(sl * yy) * cos(sm * zz)) \
          * exp(- lamb ** 2 * t * nu)                                 # lint:ok
        v = (A / (sk ** 2 + sl ** 2)) \
          * (lamb * sk * sin(sk * xx) * cos(sl * yy) * sin(sm * zz)
           - sm * sl * cos(sk * xx) * sin(sl * yy) * cos(sm * zz)) \
          * exp(- lamb ** 2 * t * nu)                                 # lint:ok
        w = A * cos(sk * xx) * cos(sl * yy) * sin(sm * zz) \
          * exp(- lamb ** 2 * t * nu)
        # del xx, yy, zz                                                # lint:ok
    elif random:
        print('random initial field')
        print('not implemented yet')
        return
    else:
        print('read initial field from file')
        print('not implemented yet')
        return

    u2uhat()  # Fourier transform initial condition

    return


def setTimeStep(tStep=tStepMax):
    global dt, intFact
    if tStepFix:
        dt = tStep
    else:
        print('Resetting time step')
        print('not implemented yet')
        return

    intFact = exp((nu * (k2xm + k2ym + k2zm) + Q) * dt)

    return


def derivatives():
    """
    Compute space derivatives from (u,v,w)hattemp
    """
    global ux, uy, uz, vx, vy, vz, wx, wy, wz
    temphat = kxm * uhattemp
    ux = np.real(ifftn(temphat))
    temphat = kym * uhattemp
    uy = np.real(ifftn(temphat))
    temphat = kzm * uhattemp
    uz = np.real(ifftn(temphat))

    temphat = kxm * vhattemp
    vx = np.real(ifftn(temphat))
    temphat = kym * vhattemp
    vy = np.real(ifftn(temphat))
    temphat = kzm * vhattemp
    vz = np.real(ifftn(temphat))

    temphat = kxm * whattemp
    wx = np.real(ifftn(temphat))
    temphat = kym * whattemp
    wy = np.real(ifftn(temphat))
    temphat = kzm * whattemp
    wz = np.real(ifftn(temphat))

    return


def nonLinear():
    """
    Compute nonlinear term for (u,v,w)hattemp
    """
    global utemp, vtemp, wtemp, nonlinuhat, nonlinvhat, nonlinwhat, temp, phat
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
    phat = -1.0 * (kxm * nonlinuhat + kym * nonlinvhat + kzm * nonlinwhat) \
           / (k2xm + k2ym + k2zm + 1.0e-13)

    nonlinuhat = - nonlinuhat - kxm * phat
    nonlinvhat = - nonlinvhat - kym * phat
    nonlinwhat = - nonlinuhat - kzm * phat

    return


def timeStep():
    global u, v, w, uhat, vhat, what, uhattemp, vhattemp, whattemp, t
    print('Time stepping')

    for n in range(1, Nt + 1):

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


def checkError():

    uErrorMax = np.max((-A / (sk ** 2 + sl ** 2)) \
              * (lamb * sl * cos(sk * xx) * sin(sl * yy) * sin(sm * zz)
                + sm * sk * sin(sk * xx) * cos(sl * yy) * cos(sm * zz)) \
              * exp(- lamb ** 2 * t * nu) - u)                        # lint:ok

    vErrorMax = np.max((A / (sk ** 2 + sl ** 2)) \
                     * (lamb * sk * sin(sk * xx) * cos(sl * yy) * sin(sm * zz)
                       - sm * sl * cos(sk * xx) * sin(sl * yy) * cos(sm * zz)) \
                     * exp(- lamb ** 2 * t * nu) - v)

    wErrorMax = np.max(A * cos(sk * xx) * cos(sl * yy) * sin(sm * zz) \
                     * exp(- lamb ** 2 * t * nu) - w)

    print('Maximum errors at the final time step are: ')
    print(uErrorMax, vErrorMax, wErrorMax)
    return