dmodelspy
---------------------------------------------------

Python conversion of Matlab dmodel routines (Battaglia et al., JVGR, 2013)

Converted
===========

- sphere: (3Dstrain surely doesn't work as there was no example to verify against)
- spheroid:
- sill: Made Sill class (v0.1.1)

ToDo
===========
- Convert sphere functions to Sphere2D and Sphere3D classes
- Convert spheroid functions to Spheroid class
- Make a Mogi class
- Add volume change to Sphere and Spheroid classes
- Add test code

Class principals
===========

Classes use the dataclass decorator

Each class should have the following functions:

- calc(x,y,z): performs all of the fundamental calculations for the given points
- calc_disp(x,y,z): returns displacements (u, v, w)
- calc_dV(x,y,z): returns volume change
- calc_tilts(x,y,z): returns tilts (East and North)
- calc_strain(x,y,z): returns strains (areal, shear1 and shear2)
- calc_all(x,y,z): returns all values: u, v, w, dV, tiltE, tiltN, areal, shear1, shear2

x, y and z can be vectors or 1d numpy arrays.  x and y should have the same
length.  z should have the same length or be a scalar, in which case the same z
will be used for each (x, y) pair)