%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Author: Alireza hadjighasem %
                        % Email: alirezah@ethz.ch     %
                        % Date:  20/7/2015            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; tic
%% Initial and final times of integration
t0 = 0;
tf = t0 + 40*(24*60*60);   % integration is equal to 40 days

Nt = 600;     % Number of intermediate times for reporting of the positions
tspan = linspace(t0,tf,Nt);  % A vector specifying the interval of integration
options = odeset('RelTol',1e-4,'AbsTol',1e-4); % ODE solver options

%% Generating a uniform grid of initial conditions
n = 200;  m = 60;
x = linspace(0,6.371*pi,n);  dx = abs(x(2)-x(1));
y = linspace(-3,3,m);        dy = abs(x(2)-x(1));
[xi,yi] = meshgrid(x,y);
%%  auxiliary distance for computing vorticity along trajectoris
rho = 0.5*dx;
%%
[xp_t,yp_t,Curlz_t] = Integrator(xi,yi,rho,tspan,options,'bickley');
%%
Curlz_avg_t = mean(Curlz_t,2);               % spatial average of vorticity
IVD = abs( Curlz_t(1,:)-Curlz_avg_t(1) );
LAVD = trapz(tspan, abs( bsxfun(@minus,Curlz_t,Curlz_avg_t) ), 1 );
%% Display LAVD OR IVD result
VMatrix = reshape(LAVD,m,n);
% VMatrix = reshape(IVD,m,n);

figure
imagesc(x,y,VMatrix);
axis equal tight;
set(gca,'ydir','normal')
%% Coherent Lagrangian vortex boundaries and vortex centers for 2D flows:
Nct = 50;                                  % Number of contour levels intended to extract
MinLength = 1.1;                           % minimal arc-length threshold
DeficiencyThresh = 0.5;                    % convexity deficiency threshold (%)
bnd = ContourExtraction(VMatrix,xi,yi,Nct,MinLength,DeficiencyThresh);
%% Vortex boundaries & centers at time t0, with the LAVD shown in the background 
figure
imagesc(x,y,VMatrix);
for kk=1:numel(bnd.xc); hold on; plot(bnd.xc{kk},bnd.yc{kk},'r','linewidth',3); end
plot(bnd.xp,bnd.yp,'or','MarkerFaceColor','r','MarkerSize',4);
axis equal tight;
set(gca,'ydir','normal')

