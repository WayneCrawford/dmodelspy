"""
Deformation forward modeling routines

==========================================================================
 USGS Software Disclaimer
 The software and related documentation were developed by the U.S.
 Geological Survey (USGS) for use by the USGS in fulfilling its mission.
 The software can be used, copied, modified, and distributed without any
 fee or cost. Use of appropriate credit is requested.

 The USGS provides no warranty, expressed or implied, as to the correctness
 of the furnished software or the suitability for any purpose. The software
 has been tested, but as with any complex software, there could be undetected
 errors. Users who find errors are requested to report them to the USGS.
 The USGS has limited resources to assist non-USGS users; however, we make
 an attempt to fix reported problems and help whenever possible.
==========================================================================
"""
from .sill import sill, sill_displ
from .sphere import sphere2D, sphere3D_displ, sphere3D_strain
from .spheroid import spheroid, spheroid_displ