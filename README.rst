Python port of Battaglia et al (2013) Matlab codes

Currently only sill/, spheroid/ and sphere/ have been converted

Also regularized the function names to {geometry}_displ (returns
displacement terms), {geometry}_strain (returns strain terms) and
{geometry} (returns all terms).  

I learned as I went along, so sill/ is probably cleaner than spheroid/, which
is cleaner than sphere/

Conversions verified using dMODELS examples, nothing more.

Still needs test cases (doctest and pytest)