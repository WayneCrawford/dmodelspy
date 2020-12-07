function dur = mctigueA8A18a(tt,nu,rho,zeta)
% sixth order free surface correction
% displacement functions (A8) and (A18) of McTigue (1987)
% Hankel transforms from Tellez et al (1997). Tables of Fourier, Laplace
% and Hankel transforms of n-dimensional generalized function. Acta 
% Applicandae Mathematicae, 48, 235-284 
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

% equation (48)
sigma = tt.^2.*exp(-tt)/(7-5*nu);
% missing Hankel transform of second part of equation (49)
tau = tt.^2.*exp(-tt)/(7-5*nu);                                


A8 = 0.5*sigma.*((1-2*nu)-tt.*zeta).*exp(-tt.*zeta).*besselj(1,tt.*rho);
A18 = 0.5*tau.*(2*(1-nu)-tt.*zeta).*exp(-tt.*zeta).*besselj(1,tt.*rho);
dur = A8 + A18;
