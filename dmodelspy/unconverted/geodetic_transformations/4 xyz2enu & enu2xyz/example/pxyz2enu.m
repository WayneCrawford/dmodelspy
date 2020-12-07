function [e,ecov]=pxyz2enu(d,dcov,origin)
%XYZ2ENU   Transforms from global cartestian to local cartesian.
%   [E,ECOV]=xyz2enu(D,DCOV,ORIGIN) transforms data vector D and
%   data covariance DCOV from a global cartesian (XYZ) coordinate
%   system to a local coordinate system aligned with the geographic
%   directions at the position ORIGIN.  D should be either 3nx1 or
%   3xn (n = number of individual vectors).  DCOV should be 3nx3n.
%   ORIGIN should be a vector of length 2 or 3.  If length 2, ORIGIN
%   is taken as a longitude, latitude pair (degrees); if length 3,
%   ORIGIN is taken as an XYZ station position. E is matrix (vector)
%   of transformed coordinates the same size as input D.  ECOV is a
%   matrix containing the transformed covariance.
%
%   E=xyz2enu(D,ORIGIN) behaves as above without a data covariance
%   matrix.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 24, 2001  Peter Cervelli        Standardized code
%   Nov 04, 2000  Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------
% USGS Software Disclaimer 
% The software and related documentation were developed by the U.S. 
% Geological Survey (USGS) for use by the USGS in fulfilling its mission. 
% The software can be used, copied, modified, and distributed without any 
% fee or cost. Use of appropriate credit is requested. 
%
% The USGS provides no warranty, expressed or implied, as to the correctness 
% of the furnished software or the suitability for any purpose. The software 
% has been tested, but as with any complex software, there could be undetected 
% errors. Users who find errors are requested to report them to the USGS. 
% The USGS has limited resources to assist non-USGS users; however, we make 
% an attempt to fix reported problems and help whenever possible. 
%==========================================================================

%Parse input arguments

    if nargin==2
        origin=dcov;
    end

    if length(origin)>2
        origin=xyz2llh(origin(:));
    end
   
    [i,j]=size(d);
    d=d(:);  

%Make transformation matrix

	Tm=pxyz2enum(origin,length(d)/3);

%Transform

    e=reshape(Tm*d,i,j);
	if nargout==2 & nargin > 2
	    ecov=Tm*dcov*Tm';
	end

