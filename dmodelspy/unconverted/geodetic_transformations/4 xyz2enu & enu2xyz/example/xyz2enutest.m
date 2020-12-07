% test script of
% function [ENU,ENUCOV]=xyz2enu(xyz,xyzcov,itrf)
% XYZ2ENU   Transforms from global cartestian to local cartesian.
%   [ENU,ENUCOV]=xyz2enu(xyz,xyzcov,itrf) transforms data vector xyz and
%   data covariance xyzcov from a global cartesian (XYZ) coordinate
%   system to a local coordinate system (ENU).
%
% =========================================================================

% =========================================================================
% TEST 1
clear all; clc; filename = 'xyz2enutest.txt';

% ITRF sites (ITRF05 frame) at epoch 2000.01.01
[SITE xyz(:,1) xyz(:,2) xyz(:,3)] = textread('augustine xyz.txt','%s %f %f %f','headerlines',1);

[ENU,LAMBDA0,PHI0,ENUCOV]=xyz2enu(xyz);

fid = fopen(filename,'w+');
fprintf(fid,'ITRF coordinates of GPS sites from Augustine volcano\n');
for i=1:size(xyz,1)
    fprintf(fid,'%s %11.3f %11.3f %11.3f \n',char(SITE(i,:)),xyz(i,:));
end
fprintf(fid,'\n');

fprintf(fid,'ENU coordinates of GPS sites from Augustine volcano\n');
fprintf(fid,'Origin (lambda0, phi0): (%11.6f %11.6f)\n',LAMBDA0,PHI0);
for i=1:size(xyz,1)
    fprintf(fid,'%s %11.3f %11.3f %11.3f \n',char(SITE(i,:)),ENU(i,:));
end
fprintf(fid,'\n');


% TEST 2
enu = ENU;
lambda0 = LAMBDA0;
phi0 = PHI0;
enucov = ENUCOV;
[XYZ,XYZCOV]=enu2xyz(enu,lambda0,phi0,enucov);

fprintf(fid,'XYZ coordinates from enu2xyz.m\n');
for i=1:size(xyz,1)
    fprintf(fid,'%s %11.3f %11.3f %11.3f \n',char(SITE(i,:)),XYZ(i,:));
end
fprintf(fid,'\n');

