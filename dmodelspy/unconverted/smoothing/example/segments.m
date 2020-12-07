function [xsi ysi xsf ysf Kodes Uss Usd deltas zst zsb] = segments(faultinp,dL,unit)
% this sscript 
%   (1) reads the fault input file fault.inp 
%                    fault.inp must be in the same folder of segment.m
%   (2) divides each fault segment in several (sub)segments of length L <= dL
%   (3) print an input file (segment.inp) with new (sub)segment parameters 
%
% INPUT
% faultinp  input file with fault parameters
% dL        a vector the same size of the number of fault segments
%               each dL(i) defines the max length of the (sub)segments for
%               the ith segment; define a single value dL to have a uniform grid
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
% SEGMENT PARAMETER
% xsi       x start
% ysi       y start
% xsf       x finish
% ysf       y finish
% zst       top
%           (positive downward and defined as depth below the reference surface)
% zsb       bottom; zb >
%           (positive downward and defined as depth below the reference surface)
% Uss/Usd   fault slip
%           strike slip fault: Uss > 0 right lateral strike slip
%           dip slip fault   : Usd > 0 reverse slip
%           tensile fault    : Uss > 0 tensile opening fault
% deltas    dip angle from horizontal reference surface (90° = vertical fault)
%           deltas can be between 0° and 90° but must be different from zero!
% *************************************************************************
%
% FAULT PARAMETER from faultinp
% xi        x start
% yi        y start
% xf        x finish
% yf        y finish
% Kode      a string that define the kind of fault: str (strike), dip or ten (tensile)
% zt        top
%           (positive downward and defined as depth below the reference surface)
% zb        bottom; zb >
%           (positive downward and defined as depth below the reference surface)
% Us/Ud     fault slip
%           strike slip fault: Us > 0 right lateral strike slip
%           dip slip fault   : Ud > 0 reverse slip
%           tensile fault    : Us > 0 tensile opening fault
% delta     dip angle from horizontal reference surface (90° = vertical fault)
%           delta can be between 0° and 90° but must be different from zero!
% *************************************************************************
%
% NOTE: the script assumes that every sub-segment has the same displacement and
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

if nargin == 2                                                              % set the default scale
    unit = 1000;
end;

% read fault input file
faultfid = fopen(faultinp,'r+');
M = textscan(faultfid,'%f %f %f %f %f %f %f %f %f %f %f %*[^\n]','headerlines',2);
xi = unit*M{2}; yi = unit*M{3}; xf = unit*M{4}; yf = unit*M{5}; Kode = M{6}; 
Us = M{7}; Ud = M{8}; delta = M{9}; zt = unit*M{10}; zb = unit*M{11};

% open segment input file
segmentfid = fopen('segments.inp','w+');
fprintf(segmentfid,'  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot\n');
fprintf(segmentfid,'xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n');


% DIVIDE FAULT IN SEGMENTS ************************************************
N = length(xi);                                                             % number of fault segments
I1 =0; I2 = 0;                                                              % initialize output vectors index

if length(dL) ~= 1 && length(dL) ~= N                                       % check that dL is either a scalar of a vector with N components
   error('dL has the wrong size'); 
end;   

for n=1:N                                                                   % loop over every segment
    
    if length(dL) == 1                                                      % assign grid step size along strike
        dS(n) = dL;
    else
        dS(n) = dL(n);
    end;
    
    segmentlength = sqrt((xi(n)-xf(n))^2+(yi(n)-yf(n))^2);                  % compute segment length
    S = ceil(segmentlength/dS(n));                                          % number of sub-segments
    
    % compute points of the grid along the segment
    xx = xi(n) + (0:1:S)*(xf(n)-xi(n))/S;
    yy = yf(n) + (yf(n)-yi(n))*(xx-xf(n))/(xf(n)-xi(n));

    XSi = xx(1:S);   YSi = yy(1:S);                                         % sub-segment initial points
    XSf = xx(2:S+1); YSf = yy(2:S+1);                                       % sub-segment final points

        % print segment.inp file
        for s=1:S
            fprintf(segmentfid,'%3d %10.1f %10.1f %10.1f %10.1f %3d %10.1f %10.1f %10.1f %10.1f %10.1f\n',...
                                n,XSi(s),YSi(s),XSf(s),YSf(s),Kode(n),Us(n),Ud(n),delta(n),zt(n),zb(n));  
        end;
        
        % assign values to output vectors
        I2 = I2 + S;                                                        % assign value to index                                                 
        I1 = I2 - S;                                                        % assign value to index
        xsi(1+I1:I2) = XSi;                                                 % x start
        ysi(1+I1:I2) = YSi;                                                 % y start
        xsf(1+I1:I2) = XSf;                                                 % x finish
        ysf(1+I1:I2) = YSf;                                                 % y finish
        Kodes(1+I1:I2) = Kode(n)*ones(size(1,S));                           % code for fault type (from Coulomb 3.1)
        Uss(1+I1:I2) = Us(n)*ones(size(1,S));                               % strike slip  
        Usd(1+I1:I2) = Ud(n)*ones(size(1,S));                               % dip slip  
        deltas(1+I1:I2) = delta(n)*ones(size(1,S));                         % dip angle  
        zst(1+I1:I2) = zt(n)*ones(size(1,S));                               % fault top  
        zsb(1+I1:I2) = zb(n)*ones(size(1,S));                               % fault bottom  
        
end;
% *************************************************************************

