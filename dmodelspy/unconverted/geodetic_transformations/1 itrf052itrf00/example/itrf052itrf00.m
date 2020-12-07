function XYZ = itrf052itrf00(xyz,epoch)
% itrf052itrf00(xyz,epoch) uses a 14-parameters Helmert's transformation to 
%  transform between ITRF05 and ITRF00. 
%  If epoch = 2000.01.01, the output is ITRF05 -> WGS 84 (G1150)
%  otherwise the function transforms    ITRF05 -> ITRF00
%
% xyz (m) is a N x 3 matrix (N = number of coordinate triples) with the 
%     ITRF05 xyz coordinates at time epoch
%              x(m)          y(m)        z(m)         ITRF site
%  xyz = [ 2612631.005	-3426807.050  4686757.927     STJO epoch 2010.04.23   
%           918129.289	-4346071.293  4561977.879]    ALGO epoch 2010.04.23
%
% epoch (yr) is the epoch of the xyz coordinates in yyyy.mm.dd, default 
%        value is 2000.01.01. Use the default value to transform ITRF05 
%        into ITRF00 = WGS 84 (G1150) 
%
% XYZ (m) is a N x 3 matrix (N = number of coordinate triples) with the 
%     ITRF00 XYZ coordinates at time epoch
%              X(m)          Y(m)        Z(m)         ITRF site
%  XYZ = [ 2612631.006  -3426807.052  4686757.923     STJO epoch 2010.04.23
%           918129.289  -4346071.295  4561977.875]    ALGO epoch 2010.04.23
%
% USAGE
%
%  xyz = [ 2612631.005	-3426807.050  4686757.927    
%           918129.289	-4346071.293  4561977.879];
%  epoch = '2010.04.23';
%
%  XYZ = itrf052itrf00(xyz,epoch);
%  disp(sprintf('%8.3f %8.3f %8.3f\n',XYZ'));
%  2612631.006 -3426807.054 4686757.908
%   918129.288 -4346071.298 4561977.860  
%
% =========================================================================
% Validated against HTDP (www.ngs.noaa.gov/TOOLS/Htdp/Htdp.shtml)
% 
% =========================================================================
% DOCUMENTATION
% The 14-parameters Helmert's transformation is (Dawson and Steed, 2004)
%
% [X00]    [dx+drx*(t-t0)]                       [x05] 
% [Y00] =  [dy+dry*(t-t0)] + (1+sc+src*(t-t0))*R*[y05]
% [Z00]    [dz+drz*(t-t0)]                       [z05]
% 
% where the rotation matrix R is the sum of a steady (Rs) and a time
% dependent (Rt) component
%
%      [ 1    Rz  -Ry]             [ 0           Rrz*(t-t0)  -Rry*(t-t0)]
% Rs = [-Rz   1    Rx]        Rt = [-Rrz*(t-t0)   0           Rrx*(t-t0)]             
%      [ Ry  -Rx   1 ]             [ Rry*(t-t0) -Rrx*(t-t0)   0         ]
%
% The transformation parameters between ITRF05 and ITRF00 are 
% (see itrf.ensg.ign.fr/ITRF_solutions/2005/tp_05-00.php)
%
% epoch     dx      dy      dz      Rx      Ry      Rz      sc
% (years)   (m)     (m)     (m)     (mas)   (mas)   (mas)   (ppb)
%  2000.0   0.0001  -0.0008 -0.0058 0.000   0.000   0.000   0.40
%   /year  -0.0002   0.0001 -0.0018 0.000   0.000   0.000   0.08
%
%
% REFERENCES
% Dawson, J. and Steed, J, International Terrestrial Reference Frame (ITRF) 
% to GDA94 Coordinate Transformations, Geoscience Australia, Geodesy's 
% Technical Paper 3795, 2004 (available on-line at 
% www.ga.gov.au/geodesy/reports/datum_transformations/)
%==========================================================================
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%   Apr 23, 2010  M. Battaglia          Original code
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

%Check input arguments
    % epoch
    if nargin==1
        t = 2000.0;
    else
        if (strcmp(epoch(5),'.')==0 || strcmp(epoch(8),'.')==0 || length(epoch)~=10) 
            error('Input t **MUST** be a string of format ''''yyyy.mm.dd''');
        end
        t1 = strcat(epoch(1:4),'.01.01');
        t2 = strcat(epoch(1:4),'.12.31');
        num = datenum(epoch,'yyyy.mm.dd')-datenum(t1,'yyyy.mm.dd');
        den = datenum(t2,'yyyy.mm.dd')-datenum(t1,'yyyy.mm.dd')+1;
        yrfrct = num/den;
        t = str2double(epoch(1:4))+yrfrct;
    end;

  % coordinates  
  if size(xyz,2)~=3
       error('Input xyz **MUST** be N x 3 matrix')
    end;

% Tranformation parameters
t0 = 2000.0;                                                % initial epoch

ds = [ 0.0001 -0.0008  -0.0058];                            % shift (m),  ds = [dx dy dz]
dr = [-0.0002  0.0001  -0.0018];                            % shift rate (m/yr), dr = [drx dry drz]
d  = [ds(1)+dr(1)*(t-t0) ds(2)+dr(2)*(t-t0) ds(3)+dr(3)*(t-t0)];

Rp = 1E-3*[0.000 0.000 0.000];                              % Rotation (as), Rp = [Rx Ry Rz]
Rrp= 1E-3*[0.000 0.000 0.000];                              % Rotation rate (as/yr), Rpr = [Rrx Rry Rrz]
Rscale = pi/(180*3600);                                     % arc-second (as) to radians (rad) 
Rp = Rscale*Rp;  Rx = Rp(1); Ry = Rp(2); Rz = Rp(3);        % scaling from as to rad 
Rrp= Rscale*Rrp; Rrx = Rrp(1); Rry = Rrp(2); Rrz = Rrp(3);  % scaling from as to rad
Rs = [ 1    Rz  -Ry 
     -Rz    1    Rx 
      Ry   -Rx   1 ];                                       % Rotation (steady)
Rt= [0              Rrz*(t-t0) -Rry*(t-t0)
     -Rrz*(t-t0)    0           Rrx*(t-t0)
      Rry*(t-t0)   -Rrx*(t-t0)  0];                         % Rotation (time-dependent)
R = Rs + Rt;                                                % Rotation matrix (rad)
  
ssc = 0.40; src = 0.08;                                     % scaling and rate (ppb)
sc = 1E-9*(ssc + src*(t-t0));                               % scaling

% Helmert's 14-parameter transformation from eq (2) of Dawson and Steed (2004)
XYZ =  ones(size(xyz,1),1)*d + ((1+sc)*R*xyz')';