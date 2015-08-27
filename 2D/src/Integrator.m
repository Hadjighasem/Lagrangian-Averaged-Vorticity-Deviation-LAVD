%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% function [xp_t,yp_t,Curlz_t] = Integrator(xi,yi,rho,tspan,options)

% Input arguments:
%   xi: x-component of a two-dimensional grid of initial conditions
%   yi: y-component of a two-dimensional grid of initial conditions
%   rho: auxiliary distance for computing vorticity along trajectoris
%   tspan: time span for advecting particles 
%   options: options structure for ordinary differential equation solvers
%   example: specify the name of the example. 'Bickley' OR 'Ocean'

% Output arguments:
%   xp_t: x-component of Lagrangian trajectories - size: [#times,#particles]
%   yp_t: y-component of Lagrangian trajectories - size: [#times,#particles]
%   Curlz_t: z-component of vorticity - size: [#times,#particles]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [xp_t,yp_t,Curlz_t] = Integrator(xi,yi,rho,tspan,options,example)
tic
% Np = numel(xi);
% Nt = numel(tspan);
%% Defining the velocity field
if strcmpi(example,'Bickley')
    
    %% Bickley jet parameters:
    U = 62.66e-6;
    L = 1770e-3;
    c2 = 0.205*U;
    c3 = 0.461*U;
    eps1 = 0.0075;
    eps2 = 0.15;
    eps3 = 0.3;
    r0 = 6371e-3;
    k1 = 2/r0; 
    k2 = 4/r0; 
    k3 = 6/r0;
    c1 = c3+((sqrt(5)-1)/2)*(k2/k1)*(c2-c3);

    u = @(x,y,t) U*sech(y/L).^2+(2*eps1*U*cos(k1*(x-c1*t))+...
                  2*eps2*U*cos(k2*(x-c2*t))+...
                  2*eps3*U*cos(k3*(x-c3*t))).*tanh(y/L).*sech(y/L).^2;

    v = @(x,y,t) -(eps1*k1*U*L*sin(k1*(x-c1*t))+...
                   eps2*k2*U*L*sin(k2*(x-c2*t))+...
                   eps3*k3*U*L*sin(k3*(x-c3*t))).*sech(y/L).^2;
elseif strcmpi(example,'Ocean')
    %% SSH data set:
   load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
   u_interp = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'linear','none');
   v_interp = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'linear','none');
   
   u = @(x,y,t) u_interp(x,y,t);
   v = @(x,y,t) v_interp(x,y,t);
else
    error('... The specified example is not valid.');
end
        
%% Computing Lagrangian trajectories:
[~,F] = ode45(@ODEfun,tspan,[xi(:);yi(:)],options,u,v);
xp_t = F(:,1:end/2);
yp_t = F(:,end/2+1:end);
%% Computing vorticity along trajectories:
Curlz_t = VorticityAlongTrajectory(xp_t,yp_t,u,v,tspan,rho);
toc
end


