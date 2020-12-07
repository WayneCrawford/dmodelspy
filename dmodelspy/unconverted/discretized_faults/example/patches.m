function patches(faultinp,dL,dW,unit)
% this script 
%   (1) reads the fault input file fault.inp 
%   (2) divides the fault segments in several patches of length <= dL and width <= dW
%   (3) print an input file (patches.inp) with the patches parameters 
%
% INPUT
% faultinp  input file with fault parameters
% dL        max length of patch 
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
% zpt       top
%           (positive downward and defined as depth below the reference surface)
% zpb       bottom; zb >
%           (positive downward and defined as depth below the reference surface)
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
segments(faultinp,dL,unit);                                                 % create segments file

% read segment input file
segmentfid = fopen('segments.inp','r+');
M = textscan(segmentfid,'%f %f %f %f %f %f %f %f %f %f %f %*[^\n]','headerlines',2);
xsi = M{2}; ysi = M{3}; xsf = M{4}; ysf = M{5}; fault = M{6}; 
Us = M{7}; Ud = M{8}; deltas = M{9}; zst = M{10}; zsb = M{11};

% open patch input file
patchfid = fopen('patches.inp','w+');
fprintf(patchfid,'  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot\n');
fprintf(patchfid,'xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n');


% DIVIDE FAULT SEGMENTS IN PATCHES ****************************************
    N = length(xsi);                                                        % number of fault segments
    delta = pi*deltas/180;                                                  % dip angle in radians
    phi = atan((xsf-xsi)./(ysf-ysi));                                       % strike angle in radians
 
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

    xpi = xxi(1:P); ypi = yyi(1:P);                                         % patch top left corner (initial point)
    xpf = xxf(1:P); ypf = yyf(1:P);                                         % patch top right corner (final point)
    zpt = zz(1:P);  zpb = zz(2:P+1);                                        % patch top and bottom 
                                        
        % print patches.inp input file
        for s=1:P
            fprintf(patchfid,'%3d %10.1f %10.1f %10.1f %10.1f %3d %10.1f %10.1f %10.1f %10.1f %10.1f\n',...
                                n,xpi(s),ypi(s),xpf(s),ypf(s),fault(n),Us(n),Ud(n),deltas(n),zpt(s),zpb(s));  
        end;
end;
% *************************************************************************
