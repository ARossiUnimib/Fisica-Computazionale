import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# A list of colors to distinguish the roots.
colors = ['b', 'r', 'g', 'y']

TOL = 1.e-8

def newton(z0, f, fprime, MAX_IT=1000):
    """The Newton-Raphson method applied to f(z).

    Returns the root found, starting with an initial guess, z0, or False
    if no convergence to tolerance TOL was reached within MAX_IT iterations.

    """

    z = z0
    for _ in range(MAX_IT):
        dz = f(z)/fprime(z)
        if abs(dz) < TOL:
            return z
        z -= dz
    return False

def plot_newton_fractal(f, fprime, n=200, domain=(2.352836, 2.352838, -0.0000005, 0.0000005)):
    """Plot a Newton Fractal by finding the roots of f(z).

    The domain used for the fractal image is the region of the complex plane
    (xmin, xmax, ymin, ymax) where z = x + iy, discretized into n values along
    each axis.

    """

    roots = []
    m = np.zeros((n, n))

    def get_root_index(roots, r):
        """Get the index of r in the list roots.

        If r is not in roots, append it to the list.

        """

        try:
            return np.where(np.isclose(roots, r, atol=TOL))[0][0]
        except IndexError:
            roots.append(r)
            return len(roots) - 1

    xmin, xmax, ymin, ymax = domain
    for ix, x in enumerate(np.linspace(xmin, xmax, n)):
        for iy, y in enumerate(np.linspace(ymin, ymax, n)):
            z0 = x + y*1j
            r = newton(z0, f, fprime)
            if r is not False:
                ir = get_root_index(roots, r)
                m[iy, ix] = ir
    nroots = len(roots)
    if nroots > len(colors):
        # Use a "continuous" colormap if there are too many roots.
        cmap = 'hsv'
    else:
        # Use a list of colors for the colormap: one for each root.
        cmap = ListedColormap(colors[:nroots])

    # Adjust extent to match the domain
    plt.figure(figsize=(8, 8))
    plt.imshow(m, cmap=cmap, origin='lower',
               extent=(xmin, xmax, ymin, ymax))  # Match the domain


    plt.axis('on')  # Turn on the axes
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')

    # Customize ticks for the axes
    specific_xticks = [2.352837350, 2.352836327, 2.352836323]
    plt.xticks(specific_xticks, rotation=45)  # Rotate ticks for readability
    plt.gca().set_xticks(specific_xticks)

    plt.savefig("newton_fractal.png")

f = lambda z: z**3 - 2*z**2 - 11*z + 12
fprime = lambda z: 3*z**2 - 4*z - 11

plot_newton_fractal(f, fprime, n=500)
