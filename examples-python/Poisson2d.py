import numpy as np

import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.ticker import LinearLocator

from jaxtyping import Float, jaxtyped
from beartype import beartype as typechecker

@jaxtyped(typechecker=typechecker)
def poisson_2d(in_data: Float[np.ndarray, "N N"], 
            N=256, 
            L=1,
            sig=0.1) -> Float[np.ndarray, "N N"]:
    
    # Translated from Matlab 
    # Originally from University of Washington, Department of Atmospheric Sciences
    # https://atmos.washington.edu/~breth/classes/AM585/lect/FS_2DPoisson.pdf
    # https://atmos.washington.edu/~breth/classes/AM585/matlab/lect/pois2Dper.m

    k = (2*np.pi/L)*(np.concat((np.arange(int(N/2)), np.arange(int(-N/2), 0))))

    kx, ky = np.meshgrid(k, k)

    delsq = -(kx**2 + ky**2)
    delsq[0, 0] = 1

    f_hat = np.fft.fft2(in_data)
    u = np.fft.ifft2(f_hat / delsq).real
    u -= u[0, 0]

    return u

def construct_f(xx, yy, N=256, L=1, sig=0.1):

    rsq = (xx-0.5*L)**2 + (yy-0.5*L)**2
    sigsq = sig**2

    return np.exp(-rsq / (2 * sigsq)) * (rsq - 2*sigsq)/(sigsq**2)

def plot_surf(f, xx, yy):

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    surf = ax.plot_surface(xx, yy, f, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


