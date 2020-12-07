function [LAMBDA0,ZONE] = lambda0(longitude)
% define the UTM central meridian LAMBDA0 and ZONE
%  Horizontal datum: WGS 84
%  Sign convention: positive E logitude; positive N latitude
%
% =========================================================================
% Validated against UTM grid map
% 
% =========================================================================
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%   May 07, 2010  M. Battaglia          Original code
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

for zone = 1:60
    lmin = -180+6*(zone-1); lmax = -180+6*zone;
    if length(longitude)==1 && longitude < lmax && longitude > lmin
        LAMBDA0 = -177+6*(zone-1);
        ZONE = zone;
    elseif (max(longitude) < lmax && min(longitude) > lmin)
        LAMBDA0 = -177+6*(zone-1);
        ZONE = zone;
    end
end    