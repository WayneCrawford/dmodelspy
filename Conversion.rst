Converting from Matlab to Python
=================================

Ran smop on m-file
----------------------
    
Removed/replaced smop-specific functions/lines
----------------------

smop always replaces `a * b` with `dot(a, b)` and `a .* b` with `multiply(a, b)`
When dealing with vectors or matrices, the matlab -> numpy equivalences are:

- Matlab ``*`` = numpy ``@`` # For matrices
- Matlab ``*`` = numpy ``*`` # If one of the values is a scalar
- Matlab ``a * b`` = numpy ``dot(a, b)`` # If it could be one or the other
- Matlab ``.*`` = numpy ``*`` 

Replaced:
##########

- ``dot(a,b)`` with ``a @ b`` (if both are always matrices)
- ``dot(a,b)`` with ``a * b`` (if a is a scalar)
- ``multiply(a,b)`` with ``a * b``
- ``fclose(all)`` with ``plt.show()``
    
Forced arrays/matrices to have known dimensions
----------------------
- for the "benchmark locations" x, y, z, start the function with
  `x = array(x, ndmin=1), `y = array(y, ndmin=1)`, `z = array(z, ndmin=1)`
