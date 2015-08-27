%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% function Curlz_t = VorticityTrajectory(xp_t,yp_t,u,v,tspan,rho)

% Input arguments:
%   xp_t: x-component of Lagrangian trajectories - size: [#times,#particles]
%   yp_t: y-component of Lagrangian trajectories - size: [#times,#particles]
%   u: an anonymous function computing the x-component of velocity at each 
%      given point in space and time ---> u(x,y,t)
%   v: an anonymous function computing the y-component of velocity at each 
%      given point in space and time ---> v(x,y,t)
%   tspan: time span used for advecting particles 
%   rho: auxiliary distance for computing vorticity along trajectoris

% Output arguments:
%   Curlz_t: z-component of vorticity - size: [#times,#particles]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function Curlz_t = VorticityAlongTrajectory(xp_t,yp_t,u,v,tspan,rho)

[Nt,Np] = size(xp_t);

Curlz_t = zeros(Nt,Np);
for kk=1:Nt
    t = tspan(kk);
    %% Computing velocity field on auxiliary grids
    Nrad = 4;
    ui = zeros(Np,Nrad);
    vi = zeros(Np,Nrad);
    for k=1:Nrad
        xp_aug = xp_t(kk,:) + rho*cos((k-1)*pi/2);
        yp_aug = yp_t(kk,:) + rho*sin((k-1)*pi/2);
    
        ui(:,k) = u(xp_aug,yp_aug,t*ones(1,Np));
        vi(:,k) = v(xp_aug,yp_aug,t*ones(1,Np));
    end
    %% Velocity gradients
%     gradx_u = ( ui(:,1)-ui(:,3) )/(2*rho);
    grady_u = ( ui(:,2)-ui(:,4) )/(2*rho);

    gradx_v = ( vi(:,1)-vi(:,3) )/(2*rho);
%     grady_v = ( vi(:,2)-vi(:,4) )/(2*rho);
    %%  Vorticity along particle trajectories
    Curlz_t(kk,:) = gradx_v-grady_u;
end
