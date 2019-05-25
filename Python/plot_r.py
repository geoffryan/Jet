import sys
import numpy as np
import matplotlib.pyplot as plt
import jetpy as jet
import grbpy as grb


def plot(filename, tS, rS, uS, thS):

    name = ".".join(filename.split("/")[-1].split("_")[-1].split('.')[:-1])

    print("Loading " + filename)
    t, r, th, prim, dV, _ = jet.loadCheckpoint(filename)

    rho = prim[:, 0]
    P = prim[:, 1]
    ur = prim[:, 2]
    ut = prim[:, 3]
    q = prim[:, 4]

    R = np.linspace(r.min(), r.max(), 4000)

    i = np.searchsorted(tS, t)

    gMax = np.sqrt(1+uS[i]*uS[i])
    RS = rS[i]

    DR = RS / (12*gMax*gMax)
    rhoMax = 4*gMax*rho0

    RHO = np.empty(R.shape)
    PPP = np.empty(R.shape)
    UR = np.empty(R.shape)

    shell = (R <= RS) & (R >= RS-DR)
    RHO[:] = rho0
    RHO[shell] = rhoMax
    UR[:] = 0.0
    UR[shell] = uS[i]
    PPP[:] = 0.0
    PPP[shell] = (gMax-1.0) * rhoMax / 3.0

    figname = "plot_r_{0:s}.png".format(name)

    fig, ax = plt.subplots(2, 3, figsize=(12, 6))
    ax[0, 0].plot(r, rho, 'k+')
    ax[0, 0].plot(R, RHO)
    ax[0, 1].plot(r, P, 'k+')
    ax[0, 1].plot(R, PPP)
    ax[0, 2].plot(r, q, 'k+')
    ax[1, 0].plot(r, ur, 'k+')
    ax[1, 0].plot(R, UR)
    ax[1, 1].plot(r, ut, 'k+')

    fig.tight_layout()
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close()

    return t, ur.max()


def calcSingleShellSolution(tmin, tmax, E0, g0, rho0, th0):

    td = np.power(9*E0 / (16*np.pi*rho0*g0*g0), 1.0/3.0)

    t = np.geomspace(3.0e-1*min(tmin, td), 3*tmax, 1000*np.log10(tmax/tmin))

    R0 = t[0] * np.sqrt(1-1.0/(g0*g0)) * grb.c
    u0 = np.sqrt(g0*g0-1)
    Mej = E0 / (g0-1) * grb.c * grb.c * grb.c

    R, u, th = grb.shock.shockEvolSpreadRK4(t, R0, u0, th0,
                                            Mej, rho0, 0.0, 0.0, 0.0, 0.0,
                                            0.0, 0.0, False)

    return t, R / grb.c, u, th


if __name__ == "__main__":

    N = len(sys.argv[1:])
    times = np.empty(N)

    for i, filename in enumerate(sys.argv[1:]):
        times[i] = jet.getTime(filename)
    times = np.array(times)

    tmin = times.min()
    tmax = times.max()

    E0 = 1.0
    g0 = 30.0
    rho0 = 1.0
    th0 = 0.5*np.pi

    tS, rS, uS, thS = calcSingleShellSolution(tmin, tmax, E0, g0, rho0, th0)

    t = np.empty(N)
    ur = np.empty(N)

    for i, filename in enumerate(sys.argv[1:]):
        t[i], ur[i] = plot(filename, tS, rS, uS, thS)

    fig, ax = plt.subplots(1, 1)
    ax.plot(t, ur)
    ax.plot(tS, uS)
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig('plot_u.png')
    plt.close()
