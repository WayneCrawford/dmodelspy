% test script of
% function LLH=xyz2llh(xyz,itrf)
% LLH=xyz2llh(xyz,itrf) calculates longitude (DD), latitude (DD), and
% height (m) from the ITRF cartestian coordinates xyz (m). 
% Geodetic positioning LLH is in the WGS 84 (G1150) reference frame.
% Horizontal datum: WGS 84; Vertical datum: WGS 84 ellipsoid.
% =========================================================================

% =========================================================================
% TEST 1 (ITRF05 to WGS84(G1150) = ITRF00 using HTDP  [epoch 2000 01 01]
clear all; clc; filename = 'xyz2llhtest.txt';

% ITRF sites (ITRF05 frame) at epoch 2000.01.01
[SITE xyz(:,1) xyz(:,2) xyz(:,3)] = textread('SiteITRF05ep2000.txt','%s %f %f %f','headerlines',1);
% ITRF sites (LLH WGS 84) at epoch 2000.01.01 from HTDP
[~, tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,4) tmp(:,5) tmp(:,6) tmp(:,7)] = ...
    textread('test htdp 2000 01 01.txt','%s %f %f %f %f %f %f %f','headerlines',1);

for i=1:length(tmp(:,1))
        htdp(:,1) = tmp(:,1)+sign(tmp(:,1)).*(tmp(:,2)/60+tmp(:,3)/3600);
        htdp(:,2) = tmp(:,4)+sign(tmp(:,4)).*(tmp(:,5)/60+tmp(:,6)/3600);
        htdp(:,3) = tmp(:,7);
end;    
     
LLH = xyz2llh(xyz);                                                         % tranformed coordinates (DD)

fid = fopen(filename,'w+');
fprintf(fid,'function LLH=xyz2llh(xyz)\n');
fprintf(fid,'see xyz2llhtest.m for detailed explanation\n');
fprintf(fid,'Comparison with NGS HTDP\n');
fprintf(fid,'\n');
for i=1:size(xyz,1)
    fprintf(fid,'Site: %s\n',char(SITE(i,:)));
    fprintf(fid,'ITRF05         (itrf): %8.3f %8.3f %9.4f\n',xyz(i,:));
    fprintf(fid,'WGS 84 (G1150) (htdp): %10.6f %10.6f %9.4f\n',htdp(i,:));
    fprintf(fid,'WGS 84 (G1150) (MatL): %10.6f %10.6f %9.4f\n',LLH(i,:));
    fprintf(fid,'MatL - htdp          : %10.6f %10.6f %9.4f\n',LLH(i,:)-htdp(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'\n');
for i=1:size(xyz,1)
    fprintf(fid,'%s %10.7f %10.7f\n',char(SITE(i,:)),LLH(i,1:2));
end


fclose(fid);