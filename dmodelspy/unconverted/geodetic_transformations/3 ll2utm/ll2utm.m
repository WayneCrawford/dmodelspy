function EN=ll2utm(ll,lambda0)
% EN=ll2utm(ll,lambda0) computes the E (east, meters) and N (north, meters)
%  coordinates on the UTM grid from the WGS 84 (G1150) longitude (DD) and
%  latitude (DD). According to DMA (1989), "The computations [...] are 
%  accurate to the nearest 0.001 arc second for geographic coordinates and 
%  to the nearest 0.01 meter for grid coordinates".
%  0.001 arc second = 3E-07 decimal degrees = 4.848E-9 rad
%  Horizontal datum: WGS 84
%  Sign convention: positive E logitude; positive N latitude
%
%
%  ll (DD) is a N x 2 matrix (N = number of coordinate doubles)
%         long(DD)    lat(DD)       ITRF site
%  ll =  [ -52.677750  47.595240     STJO
%          119.624983  49.322618]    DRAO
%
%  EN (m) is a N x 2 matrix (N = number of coordinate doubles) with the
%  east and north UTM grid coordinates
%          east (m)    north (m)    ITRF site
%  EN = [ 373873.147 5272678.216    STJO
%         690744.804 5466635.578]   DRAO
%
% USAGE
%
% see test file llh2utmtest.m
%
% 
% =========================================================================
% Validated against NGA UTM utility and ESRI ArcGIS
% 
% =========================================================================
% REFERENCES
% Defense Mapping Agency (1989). The Universal Grids: Universal Transverse
%   Mercator (UTM) and Universal Polar Stereographic (UPS). DMA Technical
%   Manual DMATM 8358.2.
% NIMA (2000), Department of Defense World Geodetic System 1984. Its 
%   Definition and Relationships with Local Geodetic Systems. Technical 
%   Report TR8350.2, 3rd edition, emendment 1.
%
% =========================================================================
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%   May 05, 2010  M. Battaglia          Original code
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


% Check input arguments
	if size(ll,2)~=2
        error('Input XYZ **MUST** be N x 2')
    end;

% Set latitude and longitude in radians  
lambda = pi*ll(:,1)/180;                                                    % longitude
lambda0= pi*lambda0/180;                                                    % central meridian
phi = pi*ll(:,2)/180;                                                       % latitude

% Set the constants for the WGS 84 ellipsoid (after NIMA, TR8350.2, 3.2)
    a=6378137;                                                              % semimajor axis (m)
    f=1/298.257223563;                                                      % flattening
% Ellipsoid parameters (DMA, 1989, 2-2.1)    
    b=(1-f)*a;                                                              % semiminor axis (m)
	e2=2*f - f^2;                                                           % eccentricity (squared)
	ep2=(a^2-b^2)/(b^2);                                                    % second eccentricity (squared)
    n = f/(2-f); n2 = n^2; n3 = n^3; n4 = n^4; n5 = n^5;
    v = a./sqrt(1-e2*sin(phi).^2);                                          % radius of curvature in the prime vertical
    Ap =           a*(1 - n+5*(n2-n3)/4 + 81*(n4-n5)/64);
    Bp =   (  3/2)*a*(n - n2 + 7*(n3-n4)/8 + 55*n5/64);
    Cp =   (15/16)*a*(n2 - n3 + 0.75*(n4-n5));
    Dp =   (35/48)*a*(n3 - n4 + 11*n5/16);
    Ep = (315/512)*a*(n4 - n5);
    S = Ap*phi-Bp*sin(2*phi)+Cp*sin(4*phi)-Dp*sin(6*phi)+Ep*sin(8*phi);     % meridional arc

% Universal Transversal Mercator projection parameters (DMA, 1989, 2-2.1)
dlambda = lambda - lambda0;                                                 % difference of longitude from central meridian
ko = 0.9996;                                                                % central scale factor
FN = zeros(size(phi)) + 10000000*(phi<0);                                   % false northing
FE = 500000;                                                                % false easting

% Terms used to calculate general equations (DMA, 1989, 2-2.3)    
ep4 = ep2^2; ep6 = ep2^3; ep8 = ep2^4;                                              
sinp = sin(phi); 
cosp = cos(phi); cos2p = cosp.^2; cos4p = cosp.^4; cos6p = cosp.^6; cos8p = cosp.^8;
tanp = tan(phi); tan2p = tanp.^2; tan4p = tanp.^4; tan6p = tanp.^6;

T1 = S*ko;
T2 = v.*sinp.*cosp*ko/2;
T3 = (T2.*cos2p/12).*(5-tan2p+9*ep2*cos2p+4*ep4*cos4p);
T4 = (T2.*cos4p/360).*(61-58*tan2p+tan4p+270*ep2*cos2p-330*tan2p*ep2.*cos2p+...
                       445*ep4*cos4p+324*ep6*cos6p-680*tan2p*ep4.*cos4p+...
                       88*ep8*cos8p-600*tan2p*ep6.*cos6p-192*tan2p.*ep8.*cos8p);
T5 = (T2.*cos6p/20160).*(1385-3111*tan2p+543*tan4p-tan6p);
T6 = v.*cosp*ko;
T7 = (T6.*cos2p/6).*(1-tan2p+ep2*cos2p);
T8 = (T6.*cos4p/120).*(5-18*tan2p+tan4p+14*ep2*cos2p-58*tan2p*ep2.*cos2p+13*ep4*cos4p+...
                     4*ep6*cos6p-64*tan2p*ep4.*cos4p-24*tan2p.*ep6.*cos6p);
T9 = (T6.*cos6p/5040).*(61-479*tan2p+179*tan4p-tan6p);

% Conversion of geographic coordinates to grid coordinates (DMA, 1989, 2-5)
N = FN + (T1 + T2.*dlambda.^2 + T3.*dlambda.^4 + T4.*dlambda.^6 + ...       % UTM grid northing (m)
               T5.*dlambda.^8);
E = FE + (T6.*dlambda + T7.*dlambda.^3 + T8.*dlambda.^5 + T9.*dlambda.^7);  % UTM grid easting (m)
EN = [E N];                                                                 % set output matrix 