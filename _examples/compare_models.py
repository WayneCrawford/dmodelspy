"""
Test 3D sill-like source, different depths
"""
from numpy import sqrt, arange, dot, finfo, zeros
import matplotlib.pyplot as plt

from dmodelspy import Sill, Spheroid, Sphere2D, Sphere3D, Mogi


eps = finfo('float').eps

# profile
x = 1000.0 * (arange(eps, 5, 0.05))
y = dot(2, x)
z = zeros(x.size)
r = sqrt(x**2 + y**2)

# common model parameters
x0, y0, z0 = 0, 0, 1000
P_G, a, nu = 0.001, 1000., 0.25

sill = Sill(x0, y0, z0, P_G=P_G, a=a, nu=nu)
sphere3D = Sphere3D(x0, y0, z0, P_G=P_G, a=a, nu=nu)
sphere2D = Sphere2D(x0, y0, z0, P_G=P_G, a=a, nu=nu)
mogi = Mogi(x0, y0, z0, dV=3e6)
spheroid = Spheroid(x0, y0, z0, P_G=P_G, a=a, nu=nu,
                    asrat=0.99, mu=9.6e9, theta=90, phi=0)


def plot_one(ax, x, y_sill, y_sphere2D, y_sphere3D, y_mogi, y_spheroid,
             ax_title):
    ax.plot(x, y_sill, 'r.', label='Sill')
    ax.plot(x, y_sphere2D, 'b.', label='Sphere2D')
    ax.plot(x, y_sphere3D, 'b-', label='Sphere3D')
    ax.plot(x, y_mogi, 'k.', label='Mogi')
    ax.plot(x, y_spheroid, 'g-', label='Spheroid')
    ax.set_title(ax_title)


usill, vsill, wsill = sill.calc_displ(x, y, 0) # PROBLEM! SILL DOESN'T TAKE Z ARRAY
u2D, v2D, w2D = sphere2D.calc_displ(x, y)
um, vm, wm = mogi.calc_displ(x, y)
u3D, v3D, w3D = sphere3D.calc_displ(x, y, 0) # HUGE TIME IF USE Z ARRAY
usphd, vsphd, wsphd = spheroid.calc_displ(x, y, z)

fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
xlim = 10
plot_one(axs[0], 0.001 * r, usill, u2D, u3D, um, usphd, 'U')
plot_one(axs[1], 0.001 * r, vsill, v2D, v3D, vm, vsphd, 'V')
plot_one(axs[2], 0.001 * r, wsill, w2D, w3D, wm, wsphd, 'W')

axs[0].set_ylabel('displacement, in meters')
axs[0].set_xlim([0, xlim])
axs[1].set_xlabel('radial distance, in kilometers')
axs[2].legend(frameon=False)
fig.suptitle('Comparison of different models')

fig, ax = plt.subplots(1, 1)
ind = arange(4)
p1 = ax.bar(ind, [sill.dV, sphere3D.dV, mogi.dV, spheroid.dV])
ax.set_xticks(ind)
ax.set_xticklabels(('Sill', 'Sphere3D', 'Mogi', 'Spheroid'))
ax.set_ylabel('dVolume [m^3]')
ax.set_title('Comparison of different models')

plt.show()

print("That's all folks!")
