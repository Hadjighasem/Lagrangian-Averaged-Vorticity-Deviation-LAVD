%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% function dy = ODEfun(t,y,u,v,w)

% Input arguments:
%   t: time
%   y: a vector of initial conditions - x, y, and z positions are concatenated 
%   u: an anonymous function computing the x-component of velocity at each 
%      given point in space and time ---> u(x,y,z,t)
%   v: an anonymous function computing the y-component of velocity at each 
%      given point in space and time ---> v(x,y,z,t)
%   w: an anonymous function computing the z-component of velocity at each 
%      given point in space and time ---> w(x,y,z,t)

% Output arguments:
%   dy: concatenated velocity vector
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function dy = ODEfun(t,y,u,v,w)

Np = numel(y)/3;           % number of particles
dy = zeros(3*Np,1);

dy(1:Np,1) = u( y(1:Np,1),y(Np+1:2*Np,1),y(2*Np+1:3*Np,1),t*ones(Np,1) );
dy(Np+1:2*Np,1) = v( y(1:Np,1),y(Np+1:2*Np,1),y(2*Np+1:3*Np,1),t*ones(Np,1) );
dy(2*Np+1:3*Np,1) = w( y(1:Np,1),y(Np+1:2*Np,1),y(2*Np+1:3*Np,1),t*ones(Np,1) );
              
end