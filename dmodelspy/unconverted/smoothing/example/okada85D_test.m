%function D = okada85D(xi,yi,xf,yf,zt,zb,delta,P,toggle)
% This script computes the smoothing operator D of a discretized rectangular fault. 
% All variable/parameters are in SI units, angles are in degrees
% 
% OUTPUT
% D         smooothing operator (Laplacian = D*U)
%
% SUB-DISLOCATION PARAMETERS (all column vectors with the same length except P and toggle) 
% xi        x start                                     [m]
% yi        y start                                     [m]
% xf        x finish                                    [m]
% yf        y finish                                    [m]
% zt        top                                         [m]
%           (positive downward and defined as depth below the reference surface)
% zb        bottom; zb > zt                             [m]
%           (positive downward and defined as depth below the reference surface)
% delta     dip angle from horizontal reference surface (90° = vertical fault)
%           delta can be between 0° and 90° but must be different from zero!
% P         number of rows 
% toggle    toggle [top bot] to define if the fault is breaking the surface [top = 1] or
%           a dike is opening [bot = 1]; default is [top bot] = [0 0];
%
% *************************************************************************
% Building the smoothing operator D for a non-uniform grid is computed in 2 steps
% 1) compute the grid vertical and horizontal increments
%
%
%            (i+1,j)
%               | 
%               | t2h
%         t1h   |      t3h
% (i,j-P)-----(i,j)-----------(i,j+P)
%               |
%               |
%               | t4h
%            (i-1,j)
%
%   assuming that the vertical increments (t2h,t4h) are constant (all the 
%   rows of the dicretized fault are equal), while the horizontal increments 
%   (t1h,t3h) may vary (the width of the fault columns can change)
%   See also Zwillnger (1997). 157. Finite difference formulas. In Handbook 
%   of Differential Equations (3rd edition), p. 661-669.   
%
% 2) compute D
% *************************************************************************

clear all; close all; path(pathdef);

faultinp = 'fault5.ss.inp';                                                 % read the input file with the fault geometry
dL = [2000 500 3000 1000 1500]; dW = 500;                                                      % set the grid size for length (dL) and width (dW) 
[xi yi xf yf U Ud delta zt zb P] = patches(faultinp,dL,dW);                  % call function

plotpatches(xi,yi,xf,yf,zt,zb,delta); hold on;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Depth (km)');
title(sprintf('Discretized fault (%d tiles)',length(xi))); 

toggle = [0 1];


% ASSIGN PARAMETERS *******************************************************
% values of the top and bottom toggle
if nargin == 8
    toggle = [0 0];
end;
    top = toggle(1);
    bot = toggle(2);

NP = length(xi);                                                            % number of tiles
N = NP/P;                                                                   % number of columns
if NP == 1
    error('single rectangular fault: cannot compute smoothing operator');
end;

delta = pi*delta/180;                                                       % dip angle in radians
phi = atan2((xf-xi),(yf-yi));                                               % strike angle in radians

% compute the location of the center point of each sub-dislocation (aka sub-fault, tile or patch)
xxm = xi + 0.5*(xf-xi);                                                     % needed to compute center of subdislocation
yym = yf + (yf-yi).*(xxm-xf)./(xf-xi);                                      % needed to compute center of subdislocation
    zbm = zt + 0.5*(zb-zt);                                                 % center of subdislocation (Up)
    xfm = xxm + (zbm-zt).*cot(delta).*cos(phi);                             % center of subdislocation (East)
    yfm = yym - (zbm-zt).*cot(delta).*sin(phi);                             % center of subdislocation (North)
xfm = xfm'; yfm = yfm'; zbm = zbm';                                         % transpose vectors
% *************************************************************************    

% COMPUTE INCREMENTS ******************************************************

% Vertical increments
t2h_ = sqrt((xfm(2)-xfm(1))^2+(yfm(2)-yfm(1))^2+(zbm(2)-zbm(1))^2);
t2h  = t2h_*ones(size(xfm));

% Horizontal increments
% Left edge (first column)
t3h(1:P) = sqrt((xfm(1+P)-xfm(1))^2+(yfm(1+P)-yfm(1))^2+(zbm(1+P)-zbm(1))^2);
t1h(1:P) = t3h(1:P);
% Right edge (last column)
t1h(NP+1-P:NP) = sqrt((xfm(NP-P)-xfm(NP))^2+(yfm(NP-P)-yfm(NP))^2+(zbm(NP-P)-zbm(NP))^2);
t3h(NP+1-P:NP) = t1h(NP+1-P:NP);
% Central columns
for n=2:N-1
    t1h(1+(n-1)*P:n*P) = sqrt((xfm(n*P-P)-xfm(n*P))^2+(yfm(n*P-P)-yfm(n*P))^2+(zbm(n*P-P)-zbm(n*P))^2);
    t3h(1+(n-1)*P:n*P) = sqrt((xfm(n*P+P)-xfm(n*P))^2+(yfm(n*P+P)-yfm(n*P))^2+(zbm(n*P+P)-zbm(n*P))^2);
end;
t1h = t1h'; t3h = t3h';                                                     % transpose vectors

% COMPUTE SMOOTHING OPERATOR D ********************************************

% 1) toggle=[0 0], fault does not break the surface and no dike intrusion                  
% compute central diagonal (0)
    d    = -2*(1./(t1h.*t3h) + 1./t2h.^2);
% compute diagonal (+1)
    dp1_ = 1./t2h.^2;
    dp1_((1:N)*P) = 0;
    dp1  = dp1_(1:NP-1);
% compute diagonal (-1)
    dm1_ = 1./t2h.^2;
    dm1_((1:N)*P+1) = 0;
    dm1  = dm1_(2:NP);
% compute diagonal (+P) 
    dp4_ = 2./(t3h.*(t1h+t3h));
    dp4  = dp4_(1:NP-P);
% compute diagonal (-P) 
    dm4_ = 2./(t1h.*(t1h+t3h));
    dm4  = dm4_(1+P:NP);

% 2) top =1, fault break the surface
if top ~= 0
    % central diagonal (0)
        d(1)             = 1./t2h(1);
        d((1:N-1)*P+1)   = 1./t2h((1:N-1)*P+1);
    % diagonal (+1)
        dp1(1)           = -1./t2h(1);
        dp1((1:N-1)*P+1) = -1./t2h((1:N-1)*P+1);
    % compute diagonal (-1)
        dm1((1:N-1)*P)   = 0;
    % compute diagonal (+P) 
        dp4(1)           = 0;
        dp4((1:N-2)*P+1) = 0;
    % compute diagonal (-P) 
        dm4(1)           = 0;
        dm4((1:N-2)*P+1) = 0;
end

% 3) bot =1, dike opening at the bottom
if bot ~= 0
    % central diagonal (0)
        d((1:N)*P)       = 1./t2h((1:N)*P);
    % compute diagonal (+1)
        dp1((1:N-1)*P)   = 0;
    % diagonal (-1)
        dm1(P-1)           = -1./t2h(P-1);
        dm1((1:N-1)*P+P-1) = -1./t2h((1:N-1)*P+P-1);
    % compute diagonal (+P) 
        dp4((1:N-1)*P)   = 0;
    % compute diagonal (-P) 
        dm4((1:N-1)*P)   = 0;
end

% smoothing operator
D = diag(d)+diag(dp1,1)+diag(dm1,-1)+diag(dp4,P)+diag(dm4,-P);
% *************************************************************************

% STEP 3: compute Laplacian
Ltmp = D*U';
N=length(xi)/P;
for n=1:N
    L(:,n) = Ltmp(1+(n-1)*P:n*P);
end;