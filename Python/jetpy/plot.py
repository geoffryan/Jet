import numpy as np
import matplotlib as mpl


def plotVarRTh(ax, f, grid, label, log=False):

    r_iph, t_jph, Index, Nr = grid

    Nth = len(t_jph)-1

    fmin = np.inf
    fmax = -np.inf

    for i in range(Nth):
        ia = Index[i]
        ib = Index[i] + Nr[i]
        fmin = min(fmin, f[ia+1:ib].min())
        fmax = max(fmax, f[ia+1:ib].max())

    if log:
        norm = mpl.colors.LogNorm(vmin=fmin, vmax=fmax)
    else:
        norm = mpl.colors.Normalize(vmin=fmin, vmax=fmax)

    for i in range(Nth):
        ia = Index[i]
        ib = Index[i] + Nr[i]

        rf = r_iph[ia:ib]
        thf = t_jph[i:i+2]

        X = rf[:, None] * np.sin(thf)[None, :]
        Z = rf[:, None] * np.cos(thf)[None, :]

        F = np.empty((Nr[i]-1, 1))
        F[:] = f[ia+1:ib, None]

        C = ax.pcolormesh(X, Z, F, norm=norm)

    ax.set_aspect('equal')
    fig = ax.get_figure()
    cb = fig.colorbar(C)
    cb.set_label(label)
    fig.tight_layout()
