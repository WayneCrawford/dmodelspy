Converting from Matlab to Python
=================================

Ran smop on m-file
---------------------------------
    
Ran 2to3 on resulting python file
---------------------------------
    
Removed/replaced smop-specific functions/lines
---------------------------------

smop replaces ``a * b`` with ``dot(a, b)`` and ``a .* b`` with ``multiply(a, b)``
When dealing with vectors or matrices, the matlab -> numpy equivalences are:

- Matlab ``*`` = numpy ``@`` # For matrices
- Matlab ``*`` = numpy ``*`` # If one of the values is a scalar
- Matlab ``a * b`` = numpy ``dot(a, b)`` # If it could be one or the other
- Matlab ``.*`` = numpy ``*`` 

Removed:
##########

- ``@function`` decorator,
- ``nargin`` and ``varargin`` lines
- ``*args`` and ``**kwargs`` attributes
- ``nargin`` attribute in function calls

Added:
##########

- ``return`` line with the appropriate attributes (had to look in original
  Matlab file)

Replaced:
##########

- ``dot(a,b)`` with ``a @ b`` (if both are always matrices)
- ``dot(a,b)`` with ``a * b`` (if a is a scalar)
- ``multiply(a,b)`` with ``a * b``
- ``fclose(all)`` with ``plt.show()``
- ``from smop import *`` with

  - ``from numpy import a, b, c ...``
  
    - where ``a, b, c`` are in the file and in the set of
    
      - ``multiply``, ``dot``, ``sin``, ``cos``, ``arange``, ``eye``, ``pi``, ``log``, ``exp``,
        ``max``, ``min``, ``abs``, ``zeros``, ``ones``, ``sinh``, ``cosh`` or other numpy functions
        that have a direct equivalence with a matlab function
      - ``finfo`` if the files contains ``eps`` 
      
        - replace ``eps`` with ``finfo(float).eps``
        
      - ``concatenate as concat``
      - ``arctan as atan``

  - ``from scipy.integrate import quad as quadl`` if ``quadl`` is in the file
  - ``from scipy.interpolate import CubicSpline as spline`` if ``spline`` is in the file
  - ``from scipy.special import jv as besselj`` if ``besselj`` is in the file

    - replace ``yt = spline(x, y, xt)`` with ``cs = spline(x, y)``, then ``yt = cs(xt)``
  
  
    
Forced arrays/matrices to have known dimensions
----------------------
- for the "benchmark locations" x, y, z, start the function with
  `x = array(x, ndmin=1), `y = array(y, ndmin=1)`, `z = array(z, ndmin=1)`

- at the end of the file, force output variables to have the same dimensionality
  as the input variables, e.g:
  `u = u.reshape(x.shape)`
