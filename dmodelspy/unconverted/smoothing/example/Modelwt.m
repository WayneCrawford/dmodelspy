function [Lap, Lap_inv] = Modelwt(nve, nhe, delx, dely, surf, varargin)
%
%	 [Lap, Lap_inv] = modelwt(nve, nhe, delx, dely, surf, {bot_displ});
%
%  Model weighting matrix to compute Laplacian and its inverse
%  for smoothing in slip inversions.
%
%  Input:
%		nve 	= number of vertical elements
%		nhe 	= number of horizontal elements
%		delx 	= length of elements in along strike dimension
%		dely 	= length of elements in dip dimension
%		surf 	= 1 if fault breaks free surface ~= 1 otherwise
%       bot_displ = (optional input) 0 if the slip at the base  
%                   of the dislocation should be tapered to zero
%                   and 1 if measurable offset should be allowed
%                   at the base (as in a dike opening)
%  Output:
%		Lap 	= finite difference Laplacian in two dimensions
%		Lap_inv = inverse of Lap

%--------------------------------------------------------------
% 5/25/00       previous file date
% 3/18/03       added some comments and switch for damping slip 
%               at base of fault (e.g. for application to dike 
%               propagation), JRM
% 10/25/04      corrected bug with derivative for bottom edge
%               in the case that bot_displ == 1;  Changed
%               ypartd( k, k-1:k) = [1 -1];
%               to
%               ypartd( k, k-1:k) = [-1 1];
%               Also changed scaling by size of subfaults for
%               surface-breaking subfaults when surf == 1 and 
%               bottom row subfaults when bot_displ == 1. It
%               used to scale them by 1/(dely^2), but since this
%               is a first (not second) difference, should just
%               be 1/dely.  -JRM
%--------------------------------------------------------------

	ngrid = nhe * nve;
	Lap = zeros(ngrid);  xpartd = zeros(ngrid); ypartd = zeros(ngrid);
    if nargin == 6
        bot_displ = varargin{1};
    else
        % If not specified, assume slip at base of 
        % dislocation tapers to zero.
        bot_displ = 0;
    end
	
%--------------------------------------------------------------
% x-derivative for central part of grid
% (exclude left & right edges)
%--------------------------------------------------------------
    for i  = nve+1:(nhe*nve-nve)
        xpartd(i,i-nve) = + 1.0;
        xpartd(i,i)     = - 2.0;
        xpartd(i,i+nve) = + 1.0;
    end
    
%--------------------------------------------------------------
% y-derivative for central part of grid
% (exclude top & bottom edges)
%--------------------------------------------------------------
	temp = zeros(nve);
	for i = 2:nve-1
		temp(i,i-1:i+1) = [1 -2 1];
	end
	for j = 1:nhe
		k = (j-1)*nve; 
		ypartd( k+1:k+nve, k+1:k+nve ) = temp;
	end

%--------------------------------------------------------------
% need to do top edge y-derivative
% for slip breaking the surface
%--------------------------------------------------------------
    for i  = 1:nhe
        k = (i-1)*nve+1;
        if surf == 1 
            % This condition should allow slip in the upper row of 
            % patchesto be nonzero when the fault breaks the surface 
            % (e.g., measurable surface offsets)
            ypartd( k, k:k+1) = dely * [-1 1];
        else
            % This condition should minimize slip in the upper row of 
            % patches (because these slip values get multiplied by -2
            % whereas the neighboring ones - in the second row from top -
            % only get multiplied by 1), thereby making the slip taper 
            % toward zero at the edge.
            ypartd( k, k:k+1) = [-2 1];
        end
    end
    
%--------------------------------------------------------------
% bottom edge y-derivative
%--------------------------------------------------------------
    for i  = 1:nhe
        k = i*nve;
        if bot_displ == 1
            ypartd( k, k-1:k) = dely * [-1 1];
        else
            ypartd( k, k-1:k) = [1 -2];
        end
    end
    

%--------------------------------------------------------------
% left edge x-derivative
%--------------------------------------------------------------
    for i  = 1:nve
        xpartd(i,i) 	    = - 2.0;
        xpartd(i,i+nve)     = + 1.0;
    end
    
%--------------------------------------------------------------
% right edge x-derivative
%--------------------------------------------------------------
    for i  = (nhe-1)*nve+1:nhe*nve
        xpartd(i,i) 	    = - 2.0;
        xpartd(i,i-nve)     = + 1.0;
    end
    
%--------------------------------------------------------------
% normalize by patch size and form inverse
%--------------------------------------------------------------
	Lap = xpartd/delx^2 + ypartd/dely^2;
	Lap_inv = inv(Lap);
