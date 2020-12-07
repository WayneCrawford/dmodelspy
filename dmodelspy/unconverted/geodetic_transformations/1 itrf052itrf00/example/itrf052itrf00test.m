% test script of
% function XYZ = itrf052itrf00(xyz,epoch)
% itrf052itrf00(xyz) uses a 14-parameters Helmert's transformation to get 
% the coordinates transformation between ITRF05 and ITRF00
%
% If epoch = 2000.01.01, the function transforms ITRF05 -> WGS 84 (G1150) = ITRF00
% otherwise the function transforms ITRF05 -> ITRF00

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


% TEST 1 (ITRF05 to WGS84(G1150) = ITRF00 using HTDP  [epoch 2000 01 01]
clear all; clc; filename = 'itrf052itrf00test.txt';

% ITRF sites (ITRF05 frame) at epoch 2000.01.01
[SITE xyz(:,1) xyz(:,2) xyz(:,3)] = textread('SiteITRF05ep2000.txt','%s %f %f %f','headerlines',1);
epoch = '2000.01.01';
% ITRF sites (ITRF00 frame) at epoch 2000.01.01 from HTDP
[~, htdp(:,1) htdp(:,2) htdp(:,3)] = textread('test htdp 2000 01 01.txt','%s %f %f %f','headerlines',1);
    
XYZ = itrf052itrf00(xyz);

fid = fopen(filename,'w+');
fprintf(fid,'function XYZ = itrf052itrf00(xyz)\n');
fprintf(fid,'see itrf052itrf00test.m for detailed explanation\n');
fprintf(fid,'\n');
fprintf(fid,'Comparison with NGS HTDP [epoch=: %s]\n',epoch);
for i=1:size(xyz,1)
    fprintf(fid,'SITE: %s\n',char(SITE(i,:)));
    fprintf(fid,'ITRF05 (itrf): %8.3f %8.3f %8.3f\n',xyz(i,:));
    fprintf(fid,'ITRF00 (htdp): %8.3f %8.3f %8.3f\n',htdp(i,:));
    fprintf(fid,'ITRF00 (MatL): %8.3f %8.3f %8.3f\n',XYZ(i,:));
    fprintf(fid,'MatL - htdp  : %8.4f %8.4f %8.4f\n',XYZ(i,:)-htdp(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');

% =========================================================================
% TEST 2 (ITRF05 to ITRF00 (WGS 84 ellipsoid) using HTDP [epoch 2010 04 23]
clear all; clc; filename = 'itrf052itrf00test.txt';

% ITRF sites (ITRF05 frame) at epoch 2010.05.14
[SITE xyz(:,1) xyz(:,2) xyz(:,3)] = textread('SiteITRF05ep2010.txt','%s %f %f %f','headerlines',1);
epoch = '2010.05.14';   
% ITRF sites (ITRF00 frame) at epoch 2010.05.14 from HTDP
[site htdp(:,1) htdp(:,2) htdp(:,3)] = textread('test htdp 2010 05 14.txt','%s %f %f %f','headerlines',1);
    
XYZ = itrf052itrf00(xyz,epoch);

fid = fopen(filename,'a+');
fprintf(fid,'\n');
fprintf(fid,'Comparison with NGS HTDP [epoch=: %s]\n',epoch);
for i=1:size(xyz,1)
    fprintf(fid,'SITE: %s\n',char(SITE(i,:)));
    fprintf(fid,'ITRF05 (itrf): %8.3f %8.3f %8.3f\n',xyz(i,:));
    fprintf(fid,'ITRF00 (htdp): %8.3f %8.3f %8.3f\n',htdp(i,:));
    fprintf(fid,'ITRF00 (MatL): %8.3f %8.3f %8.3f\n',XYZ(i,:));
    fprintf(fid,'MatL - htdp  : %8.4f %8.4f %8.4f\n',XYZ(i,:)-htdp(i,:));
    fprintf(fid,'\n');
end
fclose(fid);