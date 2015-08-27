%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% function [xp,yp,valp] = Find2DPeak(VMatrix,x,y,mode)

% Input arguments:
%   VMatrix: input scalar field.
%   x: x-component of grid points defining the input scalar field - vector
%   y: y-component of grid points defining the input scalar field - vector
%   mode: 'maxima' OR 'minima'

% Output arguments:
%   xp: x-coordinates of the local maxima or minima 
%   yp: y-coordinates of the local maxima or minima 
%   valp: valur of local minima or maxima
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [xp,yp,valp] = Find2DPeak(VMatrix,x,y,mode)
%% checks:
if strcmpi(mode,'maxima')
    interp_obj = griddedInterpolant({x,y}, permute(VMatrix,[2,1]),'spline','none');
elseif strcmpi(mode,'minima')
    interp_obj = griddedInterpolant({x,y}, permute(-VMatrix,[2,1]),'spline','none');
else
    error('The mode argument is not valid. Use minima or maxima as mode arguments');
end
%%
dx = abs( x(2)-x(1) );
dy = abs( y(2)-y(1) );
[xi,yi] = meshgrid(x,y);
[z1,z2] = gradient(VMatrix,dx,dy);

%% Step 1: Finding critical points 
%%% Surfaces crossings:
%%% Take the difference between the two surface heights and find the contour
%%% where that surface is zero.
zdiff = z1 - z2;
C = contours(xi,yi,zdiff,[0 0]);

%%% Extract the x- and y-locations from the contour matrix C.
xL = C(1,2:end);
yL = C(2,2:end);

%%% Interpolate on the first surface to find z-locations for the intersection
%%% line.
zL = interp2(xi,yi,z1,xL,yL);

%%% finding the intersection with the plane z3=0. %%%
ind = find( zL(1:end-1).*zL(2:end) <=0 );
Xc = xL(ind)+( xL(ind+1)-xL(ind) ).*( 0-zL(ind) )./( zL(ind+1)-zL(ind) );
Yc = yL(ind)+( yL(ind+1)-yL(ind) ).*( 0-zL(ind) )./( zL(ind+1)-zL(ind) );
Nc = numel(Xc);
%% Step 2: Determining the type of critical points:
%%% 1. detH > 0 and fxx(x0,y0) > 0 implies (x0,y0) is a local minimum;
%%% 2. detH > 0 and fxx(x0,y0) < 0 implies (x0,y0) is a local maximum;
%%% 3. detH < 0 implies (x0,y0) is a saddle point;
%%% 4. detH = 0 then the test is inconclusive. 
[z1,z2] = gradient(VMatrix,dx,dy);
[H11,H12] = gradient(z1,dx,dy);
[H21,H22] = gradient(z2,dx,dy);
trH = H11+H22;
detH = H11.*H22-H12.*H21;

l1 = 0.5*trH-sqrt((0.5*trH).^2-detH);
l2 = 0.5*trH+sqrt((0.5*trH).^2-detH);
l1_interp = griddedInterpolant({x,y}, permute(l1,[2,1]),'cubic','none');
l2_interp = griddedInterpolant({x,y}, permute(l2,[2,1]),'cubic','none');

xp = [];  yp = [];  valp = [];
for kk=1:Nc
    Xc_l1 = l1_interp(Xc(kk),Yc(kk));
    Yc_l2 = l2_interp(Xc(kk),Yc(kk));
    if (Xc_l1<0) && (Yc_l2<0) % local maxima
        xp(end+1) = Xc(kk);
        yp(end+1) = Yc(kk);
        valp(end+1) = interp_obj(Xc(kk),Yc(kk));
    end
end
