%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Author: Alireza hadjighasem %
                        % Email: alirezah@ethz.ch     %
                        % Date:  20/7/2015            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
load('Ocean_geostrophic_velocity.mat','date_str','time','Domain')
% The altimeter products used in this work are produced by SSALTO/DUACS and
% distributed by AVISO, with support from CNES (http://www.aviso.oceanobs.com/duacs). 
%% Initial and final times of integration
t0 = time(1);
tf = t0+90;

Nt = 600;     % Number of intermediate times for reporting of the positions
tspan = linspace(t0,tf,Nt);  % A vector specifying the interval of integration
options = odeset('RelTol',1e-4,'AbsTol',1e-4); % ODE solver options
memo1 = ['... Integration time is ',num2str(tf-t0),' days'];
%% Generating a uniform grid of initial conditions
n = 390;  m = 210;
x = linspace(-4,9,n);        dx = abs(x(2)-x(1));
y = linspace(-35,-28,m);     dy = abs(y(2)-y(1));
[xi,yi] = meshgrid(x,y);

%%  auxiliary distance for computing vorticity along trajectoris
rho = 0.5*dx;
%%
[xp_t,yp_t,Curlz_t] = Integrator(xi,yi,rho,tspan,options,'ocean');
%%
Curlz_avg_t = mean(Curlz_t,2);               % spatial average of vorticity
IVD = abs( Curlz_t(1,:)-Curlz_avg_t(1) );
LAVD = trapz(tspan, abs( bsxfun(@minus,Curlz_t,Curlz_avg_t) ), 1 );
%% Display LAVD OR IVD result
VMatrix = reshape(LAVD,m,n);
% VMatrix = reshape(IVD,m,n);

figure
imagesc(x,y,VMatrix);
set(gca,'ydir','normal')
axis equal tight; colorbar
%% Coherent Lagrangian vortex boundaries and vortex centers for 2D flows:
Nct = 50;                                  % Number of contour levels intended to extract
MinLength = 1.15;                           % minimal arc-length threshold
DeficiencyThresh = 1;                      % convexity deficiency threshold (%)
bnd = ContourExtraction(VMatrix,xi,yi,Nct,MinLength,DeficiencyThresh);
%% Vortex boundaries & centers at time t0, with the LAVD shown in the background 
figure
imagesc(x,y,VMatrix);
for kk=1:numel(bnd.xc); hold on; plot(bnd.xc{kk},bnd.yc{kk},'r','linewidth',3); end
plot(bnd.xp,bnd.yp,'or','MarkerFaceColor','r','MarkerSize',4);
axis equal tight;
set(gca,'ydir','normal')

