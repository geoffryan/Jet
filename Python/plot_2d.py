import sys
import matplotlib.pyplot as plt
import jetpy as jet


def makePlots(filename):

    name = ".".join(filename.split("/")[-1].split("_")[-1].split('.')[:-1])

    print("Loading " + filename)
    t, r, th, prim, dV, grid = jet.loadCheckpoint(filename)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    jet.plotVarRTh(ax, prim[:, 0], grid, r'$\rho_0$', True)
    figname = "plot2d_{0:s}_rho.png".format(name)
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    jet.plotVarRTh(ax, prim[:, 1], grid, r'$P$')
    figname = "plot2d_{0:s}_P.png".format(name)
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    jet.plotVarRTh(ax, prim[:, 2], grid, r'$u^{\hat{r}}$')
    figname = "plot2d_{0:s}_ur.png".format(name)
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    jet.plotVarRTh(ax, prim[:, 3], grid, r'$u^{\hat{\theta}}$')
    figname = "plot2d_{0:s}_ut.png".format(name)
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)


if __name__ == "__main__":

    for filename in sys.argv[1:]:
        makePlots(filename)
