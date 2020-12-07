function [xpi ypi xpf ypf Ups Upd deltap zpt zpb P] = patches(faultinp,dL,dW,unit)
% this script 
%   (1) reads the fault input file fault.inp 
%   (2) divides the fault segments in several patches of length <= dL and width <= dW
%   (3) print an input file (patches.inp) with the patches parameters 
%
% INPUT
% faultinp  input file with fault parameters
% dL        a vector the same size of the number of fault segments
%               each dL(i) defines the max length of the (sub)segments for
%               the ith segment; define a single value dL to have a uniform
%               grid
% dW        max width of patch 
% unit      1 if fault coordinates are in m
%           1000 if the fault units are in km [default]        
%           unit is used to scale the fault dimensions from km (Coulomb 3.1
%           standard) to m
%
% faultinp (see example below) use a fault format similar to that of the .inp files of
% Coulomb 3.1. Note that the dimensions may be m or km (Coulomb 3.1 uses km). 
% Kode is retained only for compatibility with Coulomb 3.1     
%
%  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot
%xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
%  1     5.8333     7.9839    18.0914    12.9839 100     1.0000     0.0000    90.0000     3.0000    15.0000    added by mouse-click
%
%
% PATCH PARAMETER (all parameters are in m)
% xpi       x start
% ypi       y start
% xpf       x finish
% ypf       y finish
% Ups/Upd   fault slip
%           strike slip fault: Ups > 0 right lateral strike slip
%           dip slip fault   : Upd > 0 reverse slip
%           tensile fault    : Ups > 0 tensile opening fault
% deltap    dip angle from horizontal reference surface (90° = vertical fault)
%           deltap can be between 0° and 90° but must be different from
%           zero!
% zpt       top
%           (positive downward and defined as depth below the reference surface)
% zpb       bottom; zb >
%           (positive downward and defined as depth below the reference surface)
% P         number of rows
% *************************************************************************
%
% SEGMENT PARAMETER from segment.inp (all parameters are in m)
% xsi       x start
% ysi       y start
% xsf       x finish
% ysf       y finish
% Kode      defines the kind of fault: 100 (strike and/or dip), 200 (tensile and strike), 300 (tensile and dip) 
% zst       top
%           (positive downward and defined as depth below the reference surface)
% zsb       bottom; zb >
%           (positive downward and defined as depth below the reference surface)
% Uss/Usd   fault slip
%           strike slip fault: Us > 0 right lateral strike slip
%           dip slip fault   : Ud > 0 reverse slip
%           tensile fault    : Us > 0 tensile opening fault
% deltas    dip angle from horizontal reference surface (90° = vertical fault)
%           delta can be between 0° and 90° but must be different from zero!
% *************************************************************************
%
% NOTE: the script assumes that every patch has the same displacement and
% dip angle of the original segment
% *************************************************************************
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


if nargin == 3                                                              % set the default scale
    unit = 1000;
end;

% build segment file
[xsi ysi xsf ysf Kodes Uss Usd deltas zst zsb] = segments(faultinp,dL,unit); % create segments file

% open patch input file
patchfid = fopen('patches.inp','w+');
fprintf(patchfid,'  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot\n');
fprintf(patchfid,'xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n');


% DIVIDE FAULT SEGMENTS IN PATCHES ****************************************
    N = length(xsi);                                                        % number of fault segments
    I2 = 0;                                                                 % initialize output vectors index
    delta = pi*deltas/180;                                                  % dip angle in radians
    phi = atan2((xsf-xsi),(ysf-ysi));                                        % strike angle in radians
 
for n=1:N                                                                   % for every segment

    x0i = xsi(n) + (zsb(n)-zst(n))*cot(delta(n))*cos(phi(n));               % patch lower left corner
    y0i = ysi(n) - (zsb(n)-zst(n))*cot(delta(n))*sin(phi(n));
    
    x0f = xsf(n) + (zsb(n)-zst(n))*cot(delta(n))*cos(phi(n));               % patch lower right corner
    y0f = ysf(n) - (zsb(n)-zst(n))*cot(delta(n))*sin(phi(n));

    segmentwidth = sqrt((xsi(n)-x0i)^2+(ysi(n)-y0i)^2+(zst(n)-zsb(n))^2);   % compute segment width
    P = ceil(segmentwidth/dW);                                              % number of patches
    
    % compute points of the grid along the left and right side of the segment
    zz = zst(n) + (0:1:P)*(zsb(n)-zst(n))/P;                                % grid along the depth
    
    xxi = x0i + (xsi(n)-x0i)*(zz-zsb(n))/(zst(n)-zsb(n));                   % grid along the left side of the segment
    yyi = y0i + (ysi(n)-y0i)*(zz-zsb(n))/(zst(n)-zsb(n));
    
    xxf = x0f + (xsf(n)-x0f)*(zz-zsb(n))/(zst(n)-zsb(n));                   % grid along the left side of the segment
    yyf = y0f + (ysf(n)-y0f)*(zz-zsb(n))/(zst(n)-zsb(n));

    XPi = xxi(1:P); YPi = yyi(1:P);                                         % patch top left corner (initial point)
    XPf = xxf(1:P); YPf = yyf(1:P);                                         % patch top right corner (final point)
    ZPt = zz(1:P);  ZPb = zz(2:P+1);                                        % patch top and bottom 
                                        
        % print patches.inp input file
        for s=1:P
            fprintf(patchfid,'%3d %10.1f %10.1f %10.1f %10.1f %3d %10.1f %10.1f %10.1f %10.1f %10.1f\n',...
                                n,XPi(s),YPi(s),XPf(s),YPf(s),Kodes(n),Uss(n),Usd(n),deltas(n),ZPt(s),ZPb(s));  
        end;
        
        % assign values to output vectors
        I2 = I2 + P;                                                        % assign value to index                                                 
        I1 = I2 - P;                                                        % assign value to index
        xpi(1+I1:I2) = XPi;                                                 % x start
        ypi(1+I1:I2) = YPi;                                                 % y start
        xpf(1+I1:I2) = XPf;                                                 % x finish
        ypf(1+I1:I2) = YPf;                                                 % y finish
        Ups(1+I1:I2) = Uss(n)*ones(size(1,P));                              % strike slip  
        Upd(1+I1:I2) = Usd(n)*ones(size(1,P));                              % dip slip  
        deltap(1+I1:I2) = deltas(n)*ones(size(1,P));                        % dip angle  
        zpt(1+I1:I2) = ZPt;                                                 % fault top  
        zpb(1+I1:I2) = ZPb;                                                 % fault bottom  

end;
% *************************************************************************
