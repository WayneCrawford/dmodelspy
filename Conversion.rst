Converting from Matlab to Python
=================================

Ran smop on m-file
----------------------
    
Removed/replaced smop-specific functions/lines
----------------------

smop always replaces `a * b` with `dot(a, b)` and `a .* b` with `multiply(a, b)`
When dealing with vectors or matrices, the matlab -> numpy equivalences are:
    - Matlab `*` = numpy `@`
    - Matlab `.*` = numpy `*`
The latter works for scalars, too, but `@` does not.  If one of the values is
always a scalar, one can replace `*` by `*`.  If it may or may not be a scalar,
one should use `dot()` (this is what smop does)

Replaced:
##########

- `dot(a,b)` with `a @ b` (if both are always matrices)
- `multiply(a,b)` with `a * b`
- `fclose(all)` with `plt.show`
    
Forced arrays/matrices to have known dimensions
----------------------
- for the "benchmark locations" x, y, z, start the function with
  `x = array(x, ndmin=1), `y = array(y, ndmin=1)`, `z = array(z, ndmin=1)`