function duz = mctigueA7A17(tt,nu,rho,zeta)
% first free surface correction
% displacement functions (A7) and (A17) of McTigue (1987)
%==========================================================================
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

R = sqrt(rho.^2+zeta^2);
sigma = 0.5*tt.*exp(-tt);
tau1 = sigma;

A7 = 0.5*sigma.*(2*(1-nu)-tt.*zeta).*exp(tt.*zeta).*besselj(0,tt.*R);
A17 = 0.5*tau1.*((1-2*nu)-tt.*zeta).*exp(tt.*zeta).*besselj(0,tt.*R);
duz = A7 + A17;
