function LLH=xyz2llh(xyz,itrf)
% LLH=xyz2llh(xyz,itrf) computes longitude (DD), latitude (DD), and
%  height (m) from the ITRF cartestian coordinates xyz (m). 
%  Geodetic positioning LLH is in the WGS 84 (G1150) reference frame.
%  Horizontal datum: WGS 84; Vertical datum: WGS 84 ellipsoid.
%  Sign convention: positive E longitude; positive N latitude
%
%  xyz (m) is a N x 3 matrix (N = number of coordinate triples)
%              x(m)          y(m)        z(m)         ITRF site
%  xyz = [ 2612631.005	-3426807.050  4686757.927     STJO    
%         -2059164.887	-3621108.406 4814432.301];    DRAO 
%
%  itrf is string containing the name of the ITRF realization
%       'itrf00' or 'ITRF00' 
%       'itrf05' or 'ITRF05' (default)
%   
%  LLH is a N x 3 matrix (N = number of coordinate triples) with the 
%     WGS 84 (G1150) LLH coordinates
%           long(DD)   lat(DD)  height(m)
%  LLH = [  -52.6778   47.5952  152.8469              STJO
%           119.6250   49.3226  541.8837]             DRAO
%
%
% USAGE
%
% xyz = [ 2612631.005	-3426807.050 4686757.927      
%         -2059164.887	-3621108.406 4814432.301];    
%  LLH = xyz2llh(xyz)
%
%  LLH =
%
%   -52.6778   47.5952  152.8469
%   -119.6250   49.3226  541.8837
% 
% =========================================================================
% Validated against HTDP (www.ngs.noaa.gov/TOOLS/Htdp/Htdp.shtml)
% 
% =========================================================================
% REFERENCES
% Hoffmann-Wellenhof, B., Lichtenegger, H. and J. Collins (1997). GPS.
%   Theory and Practice. 4th revised edition. Springer, New York, pp. 389
% NIMA (2000), Department of Defense World Geodetic System 1984. Its 
%   Definition and Relationships with Local Geodetic Systems. Technical 
%   Report TR8350.2, 3rd edition, emendment 1.
% NIMA (2002), Implementation of the World Geodetic System 1984 (WGS 84)
%   Reference Frame G1150. Addendum to NIMA TR 8350.2.
%
% =========================================================================
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%   May 20, 2010  M. Battaglia          Fixed minor bug in line 68
%   Apr 23, 2010  M. Battaglia          New documentation
%                                       Validation against HTDP
%                                       Updated code
%   Aug 20, 2001  Peter Cervelli        Standardized code
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
	if size(xyz,2)~=3
        error('Input XYZ **MUST** be N x 3')
    end;
    
    if nargin==1
        itrf = 'ITRF05';
    else     
         if (strcmpi(itrf,'ITRF00')==0) && (strcmpi(itrf,'ITRF05')==0) 
           error('itrf **MUST** be either ITRF00 or ITRF05')
         end
    end

% Transform ITRF05 coordinates into ITRF00 coordinates    
    if strcmpi(itrf,'ITRF05')==1 
        xyz = itrf052itrf00(xyz);
    end;
    
% Set the constants for the WGS 84 ellipsoid (after NIMA, TR8350.2, 3.2)
    a=6378137;                                                              % semimajor axis (m)
    f=1/298.257223563;                                                      % flattening
   
% Calculate longitude (rad), latitude (rad), and height (m)
% after Hoffmann-Wellenhof et al. (1997), equations (10.11) - (10.13)
    b=(1-f)*a;                                                              % semiminor axis (m)
	e2=2*f - f^2;                                                           % eccentricity (squared)
	ep2=(a^2-b^2)/(b^2);                                                    % second numerical eccentricity (10.13)
    p=sqrt(xyz(:,1).^2 + xyz(:,2).^2);                                      % radius of a parallel (10.3)

    LLH(:,1) = -atan2(xyz(:,2),xyz(:,1));                                   % longitude (rad) (10.11)
    
    theta = atan((xyz(:,3)*a)./(p*b));                                      % (10.12)
    dumb = (xyz(:,3)+ep2*b*sin(theta).^3)./(p-e2*a*cos(theta).^3);
    LLH(:,2)=atan(dumb);                                                    % latitude (rad) (10.11)
    
    N=a./sqrt(1-e2*sin(LLH(:,2)).^2);                                       % radius of curvature (10.2)
    LLH(:,3)=p./cos(LLH(:,2))-N;                                            % height (m) 

% Convert to decimal degrees (DD) and change sign so that E longitude is
% positive
    LLH(:,1:2)=LLH(:,1:2)*(180/pi)*[-1 0; 0 1];    

    