%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Author: Alireza hadjighasem %
                        % Email: alirezah@ethz.ch     %
                        % Date:  26/8/2015            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
load('SOSE_velocity.mat','date_str','time','Domain')
% The Southern Ocean State Estimation (SOSE) data set used in this work is
% produced by Matthew R. Mazloff, Patrick Heimbach, and Carl Wunsch (2010).
%% Initial and final times of integration
t0 = time(1);
tf = time(end);

Nt = 50;     % Number of intermediate times for reporting of the positions
tspan = linspace(t0,tf,Nt);  % A vector specifying the interval of integration
options = odeset('RelTol',1e-3,'AbsTol',1e-3); % ODE solver options
memo1 = ['... Integration time is ',num2str(tf-t0),' days'];
%% Generating a uniform grid of initial conditions
nx = 50;  ny = 60;  nz = 60;
x = linspace(11,16,nx);             dx = abs(x(2)-x(1));
y = linspace(-37,-33,ny);           dy = abs(y(2)-y(1));
z = linspace(7,2000,ny)*1e-3;       dz = abs(z(2)-z(1));
[xi,yi,zi] = meshgrid(x,y,z);

%%  auxiliary distances for computing vorticity along trajectoris
rho.x = 0.5*dx;   rho.y = 0.5*dy;   rho.z = 0.5*dz;
%%
[xp_t,yp_t,zp_t,Curlz_t] = Integrator(xi,yi,zi,rho,tspan,options);
%%
Curlz_avg_t = mean(Curlz_t,2);               % spatial average of vorticity

IVD = abs( Curlz_t(1,:)-Curlz_avg_t(1) );
LAVD = trapz(tspan, abs( bsxfun(@minus,Curlz_t,Curlz_avg_t) ), 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display LAVD OR IVD result
VMatrix = reshape(LAVD,ny,nx,nz);

%- Initial vortex boundaries at the base layer
BaseLayerIndex = 1;                     % isocontour based on value of outermost contour of BaseLayer
Nct = 25;                               % Number of contour levels intended to extract
MinLength = 2;                          % minimal arc-length threshold
DeficiencyThresh = 1e-2;                % convexity deficiency threshold (%)
bnd = ContourExtraction(VMatrix(:,:,BaseLayerIndex),x,y,Nct,MinLength,DeficiencyThresh);

%- Vortex boundaries & centers at time t0, with the LAVD shown in the background 
figure
imagesc(x,y,VMatrix(:,:,BaseLayerIndex));
for kk=1:numel(bnd.xc); hold on; plot(bnd.xc{kk},bnd.yc{kk},'r','linewidth',3); end
plot(bnd.xp,bnd.yp,'or','MarkerFaceColor','r','MarkerSize',4);
axis equal tight; set(gca,'ydir','normal','fontsize',12)

%- IsoSurface:
fv = isosurface(x,y,z*1e3,VMatrix,bnd.cval);

figure
p = patch(fv);
isonormals(x,y,z*1e3,VMatrix,p)
set(p,'FaceColor',[0,0.59,0.59],'EdgeColor','none');
camlight; lightangle(-94,12); lighting phong
set(gca,'zdir','reverse','color',[151,232,255]/255,'fontsize',12); 
xlabel('Longitude [\circ]'); ylabel('Latitude [\circ]');  zlabel('depth [m]');
axis tight; daspect([1,1,500]); view(3)

