import sys
import math
import numpy as np
import scipy.integrate as spintegrate
import matplotlib.pyplot as plt
import jetpy as jet
import grbpy as grb

c = 2.99792458e10
mp = 1.6726219e-24
me = 9.10938356e-28
qe = 4.803e-10
sigmaT = 6.65e-25
cgs2mJy = 1e26

verbose = True


def getTobsBounds(filenames, plot=False):

    N = len(filenames)
    t = np.empty(N)
    rmin = np.empty(N)
    rmax = np.empty(N)

    for i in range(N):
        t[i] = jet.getTime(filenames[i])
        rmin[i], rmax[i] = jet.getRminmax(filenames[i])

    if plot:
        tobsBack = t-rmin
        tobsFront = t-rmax
        tobsEq = t

        fig, ax = plt.subplots(1, 1)
        ax.plot(t, tobsBack)
        ax.plot(t, tobsFront)
        ax.plot(t, tobsEq)

        ax.set_xscale('log')
        ax.set_yscale('log')

    # Ignoring Counterjet
    tobsMin = t.min()
    tobsMax = (t-rmax).max()

    return tobsMin, tobsMax


def calcEmissivity(r, th, phi, nu, thV, t, rho, P, ur, ut, q,
                   p, epse, epsB, xiN):

    u = math.sqrt(ur*ur + ut*ut)
    g = math.sqrt(1 + u*u)
    br = ur / g
    bt = ut / g

    st = math.sin(th)
    ct = math.cos(th)
    cp = math.cos(phi)
    cv = math.cos(thV)
    sv = math.sin(thV)

    betaMu = (br*st+bt*ct)*sv*cp + (br*ct - bt*st)*cv
    if betaMu >= 1.0:
        print("WHOA THATS A BIG betaMu: {0:.10e}".format(betaMu))

    # Doppler Factor
    d = 1.0 / (g * (1-betaMu))

    # Trans-rel EoS, P/rho & internal energy
    Th = P/rho
    eth = rho * (1.5*Th + 2.25*Th*Th / (1 + math.sqrt(1+2.25*Th*Th)))

    B = math.sqrt(8*np.pi*epsB*eth)

    n = rho / mp

    gm = (p-2) * epse * eth / ((p-1) * xiN * me * n * c*c)
    gc = 3 * me * c * g / (4 * sigmaT * epsB * eth * t)

    nup = nu / d
    num = 3.0 * gm*gm * qe*B / (4*np.pi * me*c)
    nuc = 3.0 * gc*gc * qe*B / (4*np.pi * me*c)

    eP = 0.5*(p-1)*math.sqrt(3) * qe*qe*qe * xiN * n * B / (me * c*c)

    f = np.empty(nu.shape)

    if num < nuc:
        D = nup < num
        G = (num <= nup) & (nup < nuc)
        H = nuc <= nup
        f[D] = eP * np.power(nup[D]/num, 1.0/3.0)
        f[G] = eP * np.power(nup[G]/num, 0.5*(1-p))
        f[H] = eP * math.pow(nuc/num, 0.5*(1-p)) * np.power(nup[H]/nuc, -0.5*p)
    else:
        E = nup < nuc
        F = (nuc <= nup) & (nup < num)
        H = num <= nup
        f[E] = eP * np.power(nup[E]/nuc, 1.0/3.0)
        f[F] = eP * np.power(nup[F]/nuc, -0.5)
        f[H] = eP * math.pow(num/nuc, -0.5) * np.power(nup[H]/num, -0.5*p)

    return f


def f_dA(xp, zp, ra, rb, cta, sta, ctb, stb, ctv, stv, vverbose=False):

    x = zp*stv + xp*ctv
    z = zp*ctv - xp*stv

    if vverbose:
        print("f_dA: ", xp/rb, zp/rb, x/rb, z/rb)

    if z <= 0:
        return 0.0

    rxz2 = x*x + z*z
    if vverbose:
        print("      ", rxz2/(rb*rb))
    if rxz2 > rb*rb:
        return 0.0
    y2rb = rb*rb - rxz2

    if rxz2 > ra*ra:
        y2ra = 0.0
    else:
        y2ra = ra*ra - rxz2

    if vverbose:
        print("     y2r = ", y2ra/(rb*rb), y2rb/(rb*rb))

    if vverbose:
        print("     y2t = ", (stb*stb*z*z - ctb*ctb*x*x)/(rb*rb*ctb*ctb),
              (sta*sta*z*z - cta*cta*x*x)/(rb*rb*cta*cta))

    if (stb*z - ctb*x <= 0.0) != (stb*z + ctb*x <= 0):
        return 0.0
    elif ctb <= 0.0:
        y2tb = np.inf
    else:
        y2tb = (stb*z - ctb*x)*(stb*z + ctb*x) / (ctb*ctb)

    if (sta*z - cta*x <= 0.0) != (sta*z + cta*x <= 0):
        y2ta = 0.0
    else:
        y2ta = (sta*z - cta*x)*(sta*z + cta*x) / (cta*cta)

    if y2ra < 0 or y2ta < 0 or y2rb < 0 or y2tb < 0:
        print(" ruhroh - {0:.6e} {1:.6e} {2:.6e} {3:.6e}".format(
              y2ra, y2ta, y2rb, y2tb))
        print("          {0:.3e} {1:.3e} {2:.3e} {3:.3f} {4:.3f} {5:.3f}"
              " {6:.3f} {7:.3f} {8:.3f}".format(zp, ra, rb, cta, sta, ctb,
                                                stb, ctv, stv))
        print("          {0:.6e} {1:.6e}".format(x, z))

    y21 = max(y2ra, y2ta)
    y22 = min(y2rb, y2tb)

    if y21 > y22:
        return 0.0

    return math.sqrt(y22)-math.sqrt(y21)


def calcPhiDA(zp, ra, rb, tha, thb, thV):

    # Calculate the geometry of the intersection between the annular section
    # and the viewing plane. zp, xp are coordinates in the cartesian system
    # tilted to align with thV.  zp is the distance from the plane to the
    # origin, and xp is the distance along the plane.  When thV = 0 they
    # align with z and x.
    #
    # Assuming throughout that 0 <= tha,thb,thV <= pi/2

    cta = math.cos(tha)
    sta = math.sin(tha)
    ctb = math.cos(thb)
    stb = math.sin(thb)
    ctv = math.cos(thV)
    stv = math.sin(thV)

    zmax = rb*cta
    zmin = ra*ctb

    if stv != 0.0:
        xp_zmax = (zp*ctv - zmax) / stv
        xp_zmin = (zp*ctv - zmin) / stv
    else:
        # We know there is an intersection, so zmin < zp < zmax
        xp_zmax = -np.inf
        xp_zmin = np.inf

    xp_rbmax = math.sqrt((rb-zp)*(rb+zp))
    xp_rbmin = -xp_rbmax

    if ctb*ctv + stb*stv != 0.0:
        xp_tbmax = zp * (stb*ctv - ctb*stv) / (ctb*ctv + stb*stv)
    else:
        xp_tbmax = np.inf

    if ctb*ctv - stb*stv != 0.0:
        xp_tbmin = -zp * (stb*ctv + ctb*stv) / (ctb*ctv - stb*stv)
    else:
        xp_tbmin = -np.inf

    xps = np.array([xp_zmax, xp_zmin, xp_rbmax, xp_rbmin, xp_tbmin, xp_tbmax])
    xps0 = xps.copy()

    xps = xps[np.isfinite(xps)]
    xps.sort()
    xpc = 0.5*(xps[:-1]+xps[1:])
    dyc = np.array([f_dA(xp, zp, ra, rb, cta, sta, ctb, stb, ctv, stv)
                    for xp in xpc])

    good = np.nonzero(dyc > 0.0)[0]
    if good.shape[0] == 0:
        print("WHAA?")
        print(zp/rb, ra/rb, rb/rb, tha, thb, thV)
        print(xps0/rb)
        print(xps/rb)
        print((zp*stv + xps*ctv)/rb)
        print((zp*ctv - xps*stv)/rb)
        print((zp*stv + xpc*ctv)/rb)
        print((zp*ctv - xpc*stv)/rb)
        print(dyc / rb)
        for xp in xpc:
            f_dA(xp, zp, ra, rb, cta, sta, ctb, stb, ctv, stv, True)
        for i in range(len(xps)-1):
            print("  zoom: xp = ", xps[i]/rb, xps[i+1]/rb)
            xpc2 = np.linspace(xps[i], xps[i+1], 10)
            dyc2 = np.array([f_dA(xp, zp, ra, rb, cta, sta, ctb, stb,
                                  ctv, stv) for xp in xpc2])
            print(xpc2/rb)
            print(dyc2/rb)

    elif good.shape[0] > 1:
        print("HuhUH?")
        print(zp/rb, ra/rb, rb/rb, tha, thb, thV)
        print(xps/rb)
        print(dyc/rb)

    xp1 = xps[good[0]]
    xp2 = xps[good[0]+1]

    out = spintegrate.quad(f_dA, xp1, xp2, args=(zp, ra, rb, cta, sta,
                                                 ctb, stb, ctv, stv),
                           epsrel=1.0e-2)
    dA = out[0]

    xp = 0.5*(xp1+xp2)
    x = zp*stv + xp*ctv
    z = zp*ctv - xp*stv
    r = rb
    ct = z / r
    cp = x / (r*math.sqrt((1-ct)*(1+ct)))
    th = math.acos(ct)
    phi = math.acos(cp)

    return r, th, phi, dA


def calcEmissivityZone(zp, ra, rb, tha, thb, nu, thV, t, rho, P, ur, ut, q,
                       p, epse, epsB, xiN):

    # (ra, rb, tha, thb) define and annular section which intersects with the
    # plane oriented towards thV a distance zp from the origin.
    # Calculate the emissivity integrated over the intersecting area.

    # Need: effective phi of intersecting area,
    #       area of intersection
    #       emissivity of the grid cell (evaluated at phi)

    r, th, phi, dA = calcPhiDA(zp, ra, rb, tha, thb, thV)

    f = calcEmissivity(r, th, phi, nu, thV, t, rho, P, ur, ur, q,
                       p, epse, epsB, xiN)

    f *= 2*dA

    return f


def calcEmissivityLayer(zp, nu, thV, t, r, th, rho, P, ur, ut, q,
                        riph, tjph, Index, Nr,
                        p, epse, epsB, xiN):

    # Ok, so this function integrates over a plane intersecting the
    # Jet grid.  The plane has normal (th, phi) = (thV, 0), and is a distance
    # zp from the origin. zp is the vertical coordinate in a coordinate
    # system oriented around (thV, 0).  Since each cell in Jet is really
    # an annulus in longitude, calculating particular intersection is a
    # little tricky.

    f = np.zeros(nu.shape)

    # Loop over all rays
    for j in range(len(tjph)-1):
        # Indices into global arrays for this ray. Ignoring first (ghost)
        # cell.
        ia = Index[j]+1
        ib = Index[j]+Nr[j]

        tha = tjph[j]
        thb = tjph[j+1]
        ariph = riph[ia-1:ib]

        camv = math.cos(tha-thV)
        cbmv = math.cos(thb-thV)
        capv = math.cos(tha+thV)
        cbpv = math.cos(thb+thV)

        ctvs = np.array([camv, cbmv, capv, cbpv])

        # Calculate max and min achievable zp's on this ray
        zpas = ariph[0] * ctvs
        zpbs = ariph[-1] * ctvs
        zpMin = min(zpas.min(), zpbs.min())
        zpMax = max(zpas.max(), zpbs.max())

        # Skip if there's no possible intersection
        if zp <= zpMin or zp >= zpMax:
            continue

        arho = rho[ia:ib]
        aP = P[ia:ib]
        aur = ur[ia:ib]
        aut = ut[ia:ib]
        aq = q[ia:ib]

        cs2iso = aP/arho
        cs2cut = 1.0e-5 * c*c

        # Calculate the zp's along each ray
        zps = ariph[:, None] * ctvs[None, :]
        zpMaxs = (zps > zp).any(axis=1)
        zpMins = (zps < zp).any(axis=1)

        # Get the cells which have an intersection & have been shocked
        good = np.nonzero((zpMins[:-1] | zpMins[1:])
                          & (zpMaxs[:-1] | zpMaxs[1:])
                          & (cs2iso > cs2cut))[0]

        if len(good) == 0:
            continue

        if verbose:
            print("        Ray {0:d} of {1:d} (Nr={2:d}): {3:d}-{4:d}"
                  .format(j, len(tjph)-1, Nr[j]-1, good[0], good[-1]))

        for i in good:
            f += calcEmissivityZone(zp, ariph[i], ariph[i+1], tha, thb,
                                    nu, thV,
                                    t, arho[i], aP[i], aur[i], aut[i], aq[i],
                                    p, epse, epsB, xiN)

    return f


def calcFluxStep(tobs, nu, filename, thV, E0, n0, p, epse, epsB, xiN):

    t, r, th, prim, dV, grid = jet.loadCheckpoint(filename)
    print("Computing contribution from {0:s}.  te = {1:.3e}".format(
          filename, t))

    riph = grid[0]
    tjph = grid[1]
    Index = grid[2]
    Nr = grid[3]

    rho = prim[:, 0]
    P = prim[:, 1]
    ur = prim[:, 2]
    ut = prim[:, 3]
    q = prim[:, 4]

    rho0 = n0*mp
    tau = math.pow(E0/(mp*n0 * c**5), 1.0/3.0)

    riph *= c*tau

    t *= tau
    r *= c*tau
    rho *= rho0
    P *= rho0*c*c

    f = np.empty((len(tobs), len(nu)))

    for i in range(len(tobs)):
        zp = c * (t - tobs[i])
        if verbose:
            print("    Computing for tobs = {0:.3e}s"
                  ", z-prime = {1:.3e}cm".format(tobs[i], zp))
        f[i, :] = calcEmissivityLayer(zp, nu, thV, t, r, th, rho, P, ur, ut,
                                      q, riph, tjph, Index, Nr,
                                      p, epse, epsB, xiN)

    return f


def calcLC(tobs, nu, filenames, thV, E0, n0, p, epse, epsB, xiN, dL, z):

    N = len(filenames)
    t = np.empty(N)

    tau = math.pow(E0/(mp*n0 * c**5), 1.0/3.0)

    for i in range(N):
        t[i] = tau * jet.getTime(filenames[i])

    sortinds = np.argsort(t)
    t = t[sortinds]

    tobsz = tobs / (1+z)
    nuz = (1+z) * np.atleast_1d(nu)

    f = np.empty((len(tobsz), len(nuz), N))

    for i in range(N):
        f[:, :, i] = calcFluxStep(tobsz, nuz, filenames[sortinds[i]],
                                  thV, E0, n0, p, epse, epsB, xiN)

    zi = c*(t[None, :] - tobsz[:, None])

    dz = zi[:, 1:] - zi[:, :-1]
    fs = 0.5*(f[:, :, :-1] + f[:, :, 1:])

    Fnu = (fs[:, :, :] * dz[:, None, :]).sum(axis=2)
    Fnu *= (1 + z) * cgs2mJy / (4 * np.pi * dL*dL)

    if len(nuz) == 1:
        Fnu = Fnu[:, 0]

    return Fnu


if __name__ == "__main__":

    filenames = sys.argv[1:]

    N = len(filenames)

    tobsMin, tobsMax = getTobsBounds(filenames)
    print(tobsMin, tobsMax)

    tobs = np.geomspace(tobsMin, tobsMax, 20)
    nu = np.array([1.0e9, 1.0e14, 1.0e19])

    thV = 0.0
    E0 = 1.0e53
    n0 = 1.0
    p = 2.2
    epse = 0.1
    epsB = 0.01
    xiN = 1.0
    dL = 1.0e28
    z = 0.0

    tau = math.pow(E0/(mp*n0 * c**5), 1.0/3.0)
    tobs *= tau

    Fnu = calcLC(tobs, nu, filenames, thV, E0, n0, p, epse, epsB, xiN, dL, z)

    f = open("lc_flux.txt", "w")
    header = "#"
    header += " thV={0:.3e}".format(thV)
    header += " E0={0:.3e}".format(E0)
    header += " n0={0:.3e}".format(n0)
    header += " p={0:.3e}".format(p)
    header += " epse={0:.3e}".format(epse)
    header += " epsB={0:.3e}".format(epsB)
    header += " xiN={0:.3e}".format(xiN)
    header += " dL={0:.3e}".format(dL)
    header += " z={0:.3e}".format(z)
    f.write(header + "\n")
    f.write("t_obs(s) nu(Hz) F_nu(mJy)\n")

    for i in range(len(nu)):
        for j in range(len(tobs)):
            line = "{0:.3e} {1:.3e} {2:.3e}\n".format(tobs[j], nu[i],
                                                      Fnu[j, i])
            f.write(line)
    f.close()

    atobs = np.geomspace(0.001*tobs[0], tobs[-1], 300)
    anu = np.empty(atobs.shape)
    thC = 0.2
    thW = thC
    Y = np.array([thV, E0, thC, thW, 0.0, 0.0, 0.0, 0.0, n0, p, epse, epsB,
                  xiN, dL, 10])
    anu[:] = nu[0]
    FnuAPr = grb.fluxDensity(atobs, anu, -1, 0, *Y)
    anu[:] = nu[1]
    FnuAPo = grb.fluxDensity(atobs, anu, -1, 0, *Y)
    anu[:] = nu[2]
    FnuAPx = grb.fluxDensity(atobs, anu, -1, 0, *Y)

    fig, ax = plt.subplots(1, 1)
    ax.plot(tobs, Fnu[:, 0], ls='-', color='C0')
    ax.plot(tobs, Fnu[:, 1], ls='-', color='C1')
    ax.plot(tobs, Fnu[:, 2], ls='-', color='C2')
    ax.plot(atobs, FnuAPr, ls='--', color='C0')
    ax.plot(atobs, FnuAPo, ls='--', color='C1')
    ax.plot(atobs, FnuAPx, ls='--', color='C2')

    Y[-1] = 30
    anu[:] = nu[0]
    FnuAPr = grb.fluxDensity(atobs, anu, -1, 0, *Y)
    anu[:] = nu[1]
    FnuAPo = grb.fluxDensity(atobs, anu, -1, 0, *Y)
    anu[:] = nu[2]
    FnuAPx = grb.fluxDensity(atobs, anu, -1, 0, *Y)
    ax.plot(atobs, FnuAPr, ls=':', color='C0')
    ax.plot(atobs, FnuAPo, ls=':', color='C1')
    ax.plot(atobs, FnuAPx, ls=':', color='C2')

    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()
