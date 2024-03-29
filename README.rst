dmodelspy
---------------------------------------------------

Python conversion of Matlab dmodels routines (Battaglia et al., JVGR, 2013)

Converted
================

- sill: Made Sill class (v0.1.1)
- sphere: Made Sphere2D and Sphere3D classes (v0.1.2)
- spheroid: Make Spheroid class (v0.1.3)

Added
================

- Mogi class (v0.1.3)

Class principals
================

Classes use the dataclass decorator

Each class should have the dependant property `dV` and the functions:

- calc_displ(x,y,z): returns displacements (u, v, w)
- calc_tilts(x,y,z): returns tilts (East and North)
- calc_strain(x,y,z): returns strains (areal, shear1 and shear2)
- calc_all(x,y,z): returns all values: u, v, w, tiltE, tiltN, areal, shear1, shear2

x, y and z can be vectors or 1d numpy arrays.  x and y should have the same
length.  z should have the same length or be a scalar, in which case the same z
will be used for each (x, y) pair)

Exceptions
**********

- Mogi and Sphere2D classes only use (x, y) on input
- Mogi class has dV as input