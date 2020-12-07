function [XYZ,XYZCOV]=enu2xyz(enu,lambda0,phi0,enucov)
% ENU2XYZ   Transforms from global cartestian to local cartesian.
%   [XYZ,XYZCOV]=enu2xyz(enu,lambda0,phi0,enucov) transforms data vector 
%   enu and data covariance enucov from a local cartesian (ENU) coordinate
%   system to the ITRF05 coordinate system (XYZ).
%   Default: 
%   enucov = I (identity matrix)
%
%  enu is a N x 3 matrix (N = number of coordinate triples) with the 
%      local ENU coordinates
%            East(m)   North(m)  Up(m)
%  enu = [ 1304312.389 -3541654.603 5707887.617
%         -4839699.176  402049.521 -2438613.305]
%
%  enucov (m^2) is a 3N x 3N matrix (N = number of coordinate triples)
%       containing the transformed covariance. Default is 
%       enucov = I (identity matrix)
%
%  (lambda0, phi0) defines the longitude and latitude of the origin of the 
%       local coordinate system in
%        
%
%  XYZ (m) is a N x 3 matrix (N = number of coordinate triples)
%              x(m)          y(m)        z(m)         
%  XYZ = [ -2900773.239 -1440890.357 5476476.756
%          -2911728.045 -1461015.567 5465156.878]
%
%  XYZCOV (m^2) is a 3N x 3N matrix (N = number of coordinate triples) 
%      containing the covariance associated to XYZ    
%
% =========================================================================
% USAGE
% see xyz2enutest.m
%
% =========================================================================
% Validated against xyz2enu.m (see xyz2enutest.m)
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
%   May 16, 2010  M. Battaglia          Original code
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
    [n,m] = size(enu);
    N = 3*n;                                                                % dimensions of covariance matrix
    if nargin==3
        enucov = eye(3*size(xyz,1),3*size(xyz,1));
    end

% Make transformation matrix
    origin=[lambda0,phi0]*pi/180;                                           % Convert to radians and evaluate trigonometric functions
    s=sin(origin);                                                          % s(1) = sin(lambda), s(2) = sin(phi)
    c=cos(origin);                                                          % c(1) = cos(lambda), c(2) = cos(phi)

% Make the (sparse) transformation matrix [see en.wikipedia.org/wiki/Geodetic_system]
    i=[1:3:N 2:3:N 1:3:N 2:3:N 3:3:N 1:3:N 2:3:N 3:3:N]; 
    j=[1:3:N 1:3:N 2:3:N 2:3:N 2:3:N 3:3:N 3:3:N 3:3:N];
    T=repmat([-s(1) c(1) -s(2)*c(1) -s(2)*s(1) c(2) c(2)*c(1) c(2)*s(1) s(2)],n,1);
    Tm=sparse(i,j,T(:));

% Transform
    enu = enu(:);
    XYZ=reshape(Tm*enu,n,m);
	XYZCOV=Tm*enucov*Tm';  
    