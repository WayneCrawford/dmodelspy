Python port of Battaglia et al (2013) Matlab codes

Currently only ``sill``, ``spheroid`` and ``sphere`` have been converted

You can install the module in your environment by by installing `pip`, then typing:
``pip install -e dmodelspy`` at the command line, just above the dmodelspy directory

The ``_examples/`` directory has examples of running the codes

Regularized the function names to 

- ``{geometry}_displ`` (returns displacement terms)
- ``{geometry}_strain`` (returns strain terms)
- ``{geometry}`` (returns all terms).  

But not all geometries have all these funtions. Current functions are ``sphere2D``, ``sphere3D_displ``, ``sphere3D_strain``,
``spheroid``, ``spheroid_displ``, ``sill`` and ``sill_displ``.

Conversions verified using dMODELS examples, nothing more.

Still needs test cases (doctest and pytest)
