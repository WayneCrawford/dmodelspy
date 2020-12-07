function [ENU,LAMBDA0,PHI0,ENUCOV]=xyz2enu(xyz,xyzcov,itrf)
%XYZ2ENU   Transforms from global cartestian to local cartesian.
%   [ENU,ENUCOV]=xyz2enu(xyz,xyzcov,itrf) transforms data vector xyz and
%   data covariance xyzcov from a global cartesian (XYZ) coordinate
%   system to a local coordinate system (ENU).
%   The origin of the local coordinate system is set to min(xyz).
%   Default: 
%   xyzcov = I (identity matrix)
%   itrf   = ITRF05
%
%  xyz (m) is a N x 3 matrix (N = number of coordinate triples)
%              x(m)          y(m)        z(m)         
%  xyz = [ -2900773.239 -1440890.357 5476476.756
%          -2911728.045 -1461015.567 5465156.878]
%
%  itrf is string containing the name of the ITRF realization
%       'itrf00' or 'ITRF00' 
%       'itrf05' or 'ITRF05' (default)
%
%  xyzcov (m^2) is a 3N x 3N matrix (N = number of coordinate triples) 
%      containing the covariance associated to xyz 
%       default xyzcov = I (identity matrix)
%   
%  ENU is a N x 3 matrix (N = number of coordinate triples) with the 
%      local ENU coordinates
%            East(m)   North(m)  Up(m)
%  ENU = [ 1304312.389 -3541654.603 5707887.617
%         -4839699.176  402049.521 -2438613.305]
%
%  ENUCOV (m^2) is a 3N x 3N matrix (N = number of coordinate triples)
%       containing the transformed covariance.
%
% =========================================================================
% USAGE
% see xyz2enutest.m
%
% =========================================================================
% Validated against enu2xyz (see xyz2enutest.m)
% 
% =========================================================================
% REFERENCES
% Hoffmann-Wellenhof, B., Lichtenegger, H. and J. Collins (1997). GPS.
%   Theory and Practice. 4th revised edition. Springer, New York, pp. 389
% 
% =========================================================================
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%   May 16, 2010  M. Battaglia          New documentation
%                                       Validation against enu2xyz
%                                       Updated code
%   Aug 24, 2001  Peter Cervelli        Standardized code
%   Nov 04, 2000  Peter Cervelli		Original Code
%
% =========================================================================
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


% Parse input arguments
    [n,m] = size(xyz);
    N = 3*n;                                                                % dimensions of covariance matrix
    if nargin==2
        itrf = 'ITRF05';
    end
 
    if nargin==1
        itrf = 'ITRF05';
        xyzcov = eye(3*size(xyz,1),3*size(xyz,1));
    end
    
    origin =xyz2llh(min(xyz),itrf);                                         % set origin
    origin = origin(1:2);
    
% Make transformation matrix
    origin=origin*pi/180;                                                   % Convert to radians and evaluate trigonometric functions
    s=sin(origin);                                                          % s(1) = sin(lambda), s(2) = sin(phi)
    c=cos(origin);                                                          % c(1) = cos(lambda), c(2) = cos(phi)

% Make the transformation matrix [Hoffmann-Wellenhof, 1987, p. 149, (7.6), see en.wikipedia.org/wiki/Geodetic_system]
    i=[1:3:N 2:3:N 3:3:N 1:3:N 2:3:N 3:3:N 2:3:N 3:3:N];
    j=[1:3:N 1:3:N 1:3:N 2:3:N 2:3:N 2:3:N 3:3:N 3:3:N];
    T=repmat([-s(1) -s(2)*c(1) c(2)*c(1) c(1) -s(2)*s(1) c(2)*s(1) c(2) s(2)],n,1);
    Tm=sparse(i,j,T(:));
    
% Transform
    xyz = xyz(:);
    ENU=reshape(Tm*xyz,n,m);
	ENUCOV=Tm*xyzcov*Tm';  
    LAMBDA0 = origin(1)*180/pi;
    PHI0 = origin(2)*180/pi;
    