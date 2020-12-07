% test script of
% function EN=llh2utm(ll,lambda0)
% EN=llh2utm(ll,lambda0) computes the E (east, meters) and N (north, meters)
%  coordinates on the UTM grid from the WGS 84 (G1150) longitude (DD) and
%  latitude (DD). According to DMA (1989), "The computations [...] are 
%  accurate to the nearest 0.001 arc second for geographic coordinates and 
%  to the nearest 0.01 meter for grid coordinates".
%  0.001 arc second = 3E-07 decimal degrees = 4.848E-9 rad
%  Horizontal datum: WGS 84
%  Sign convention: positive E logitude; positive N latitude
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

% =========================================================================
% TEST 1 
clear all; clc; filename = 'll2utmtest.txt';

% ITRF sites (ITRF05 frame) at epoch 2000.01.01
[SITE xyz(:,1) xyz(:,2) xyz(:,3)] = textread('SiteITRF05ep2000.txt','%s %f %f %f','headerlines',1);
% ITRF sites (UTM WGS 84) at epoch 2000.01.01 from ArcGIS
[site, ~, ~, utm(:,1) utm(:,2)] = textread('test ArcGIS UTM.txt','%s %f %f %f %f','headerlines',1);

LLH = xyz2llh(xyz);                                                         % tranformed coordinates (DD)

LAMBDA0 = zeros(size(LLH(:,1)));                                            % variable pre-allocation
ZONE = zeros(size(LLH(:,1)));                                               % variable pre-allocation
EN = zeros(size(LLH(:,1:2)));                                               % variable pre-allocation

for i=1:length(LLH(:,1))
    [LAMBDA0(i),ZONE(i)] = lambda0(LLH(i,1));
    EN(i,:) = ll2utm(LLH(i,1:2),LAMBDA0(i));
end



fid = fopen(filename,'w+');
fprintf(fid,'function EN=ll2utm(ll,lambda0)\n');
fprintf(fid,'see ll2utmtest.m for detailed explanation\n');
fprintf(fid,'Comparison with ArcGIS\n');
fprintf(fid,'\n');

for i=1:size(EN,1)
    fprintf(fid,'SITE: %s %s\n',char(SITE(i,:)),char(site(i,:)));
    fprintf(fid,'WGS 84 (G1150) (MatL): %11.7f %11.7f\n',LLH(i,1:2));
    fprintf(fid,'UTM (WGS 84) (ArcGIS): %11.3f %11.3f\n',utm(i,:));
    fprintf(fid,'UTM (WGS 84) (MatLab): %11.3f %11.3f %4.1f %3i\n',EN(i,:),LAMBDA0(i),ZONE(i));
    fprintf(fid,'MatLab - ArcGIS      : %11.3f %11.3f\n',EN(i,:)-utm(i,1:2));
    fprintf(fid,'\n');

end