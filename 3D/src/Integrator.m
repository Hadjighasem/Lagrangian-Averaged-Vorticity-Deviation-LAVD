%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% function [xp_t,yp_t,zp_t,Curlz_t] = Integrator(xi,yi,zi,rho,tspan,options)

% Input arguments:
%   xi: x-component of a two-dimensional grid of initial conditions
%   yi: y-component of a two-dimensional grid of initial conditions
%   zi: z-component of a two-dimensional grid of initial conditions
%   rho: auxiliary distance for computing vorticity along trajectoris
%   tspan: time span for advecting particles 
%   options: options structure for ordinary differential equation solvers

% Output arguments:
%   xp_t: x-component of Lagrangian trajectories - size: [#times,#particles]
%   yp_t: y-component of Lagrangian trajectories - size: [#times,#particles]
%   zp_t: z-component of Lagrangian trajectories - size: [#times,#particles]
%   Curlz_t: z-component of vorticity - size: [#times,#particles]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [xp_t,yp_t,zp_t,Curlz_t] = Integrator(xi,yi,zi,rho,tspan,options)
tic
% Np = numel(xi);
% Nt = numel(tspan);
%% Defining the velocity field
%- SOSE data set:
load('SOSE_velocity.mat','lon','lat','depth','UT','VT','WT','time');
u_interp = griddedInterpolant({lon,lat,depth,time},permute(UT,[2,1,3,4]),'linear','nearest');
v_interp = griddedInterpolant({lon,lat,depth,time},permute(VT,[2,1,3,4]),'linear','nearest');
w_interp = griddedInterpolant({lon,lat,depth,time},permute(WT,[2,1,3,4]),'linear','nearest');

u = @(x,y,z,t) u_interp(x,y,z,t);
v = @(x,y,z,t) v_interp(x,y,z,t);
w = @(x,y,z,t) w_interp(x,y,z,t);
        
%% Computing Lagrangian trajectories:
[~,F] = ode45(@ODEfun,tspan,[xi(:);yi(:);zi(:)],options,u,v,w);
xp_t = F(:,1:end/3);
yp_t = F(:,end/3+1:2*end/3);
zp_t = F(:,2*end/3+1:end);
%% Computing vorticity along trajectories:
%- NOTE: The x and y components of vorticity are negligible
Curlz_t = VorticityAlongTrajectory(xp_t,yp_t,zp_t,u,v,tspan,rho);
toc
end


