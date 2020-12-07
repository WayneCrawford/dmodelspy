sphere/
--------------
converted to python, example works.

    - mctigue3Dstrain surely doesn't work, as there is no example so I
      didn't apply the same changes/verifications as I did to mctigue3Ddispl

spheroid:
--------------
partially converted to python
    - ran smop
    - ran 2to3
    - did basic smop emancipation
        - replace smop module call by numpy routines
        - remove @function decorator, nargin and varargin lines, *args and
          **kwargs attributes, added return with attributes
Need to:
    - run matlab examples and save results
    - modify python example codes, then run and compare results