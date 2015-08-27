%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% bnd = ContourExtraction(VMatrix,xi,yi,Nct,MinLength,DeficiencyThresh)

% Input arguments:
%   VMatrix: input scalar field.
%   xi: x-component of grid points defining the input scalar field - matrix
%   yi: y-component of grid points defining the input scalar field - matrix
%   Nct: Number of contour levels intended to extract
%   MinLength: minimal arc-length threshold
%   DeficiencyThresh: convexity deficiency threshold - as a percentage of
%   convex hull

% Output arguments:
%   bnd: a structure containing positions of vortex boundaries and vortex
%   centers - struct
%       bnd.xc: x-coordinates of contours defining vortex boundaries - cell
%       bnd.yc: y-coordinates of contours defining vortex boundaries - cell
%       bnd.cval: corresponding contour values on the scalar field - vector
%       bnd.xp: x-coordinates of vortex centers - vector
%       bnd.yp: y-coordinates of vortex centers - vector
%       bnd.valp: corresponding values of vortex cores on the scalare field - vector
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function bnd = ContourExtraction(VMatrix,xi,yi,Nct,MinLength,DeficiencyThresh)
x = xi(1,:);
y = yi(:,1);  
%% Step1: Extracting closed, (almost) convex contours
h1 = contourc(x,y,VMatrix,Nct);
s = getcontourlines(h1);
disp(sprintf('... %3d contours are extracted.',numel(s)));
ct_x = [];   ct_y = [];   ct_val = [];
for k=1:numel(s)
    % check if the contour is closed
    if ( s(k).x(1)==s(k).x(end) ) && ( s(k).y(1)==s(k).y(end) ) 
        Length = sum( sqrt(diff(s(k).x,1).^2+diff(s(k).y,1).^2) );
        if Length>MinLength
            % Check if the contour is convex
            if IsContourConvex(s(k).x,s(k).y,DeficiencyThresh) 
                ct_x{end+1} = s(k).x';
                ct_y{end+1} = s(k).y'; 
                ct_val{end+1} = s(k).v';
            end
        end
    end
end
%% Step2: Selecting the outermost contour for each nested family of contours
Nct_filt = numel(ct_x);
disp(sprintf('... %3d contours are closed and convex.',Nct_filt));

% selecting the first point of each contour as a query point
if Nct_filt~=0
    xq = cellfun(@(x) x(1),ct_x,'UniformoutPut',true);  
    yq = cellfun(@(x) x(1),ct_y,'UniformoutPut',true);  
end

%%
indEliminte_max = false(Nct_filt,1);
for ii=1:Nct_filt
    in = inpolygon(xq,yq,ct_x{ii},ct_y{ii});  in(ii)=0;
    indEliminte_max(in) = 1;
end
ct_x = ct_x(~indEliminte_max);
ct_y = ct_y(~indEliminte_max);
ct_val = ct_val(~indEliminte_max);

%% Step3: Finding vortex centers
%%% Step1: extracting local maxima and minima of the field
[xp,yp,valp] = Find2DPeak(VMatrix,x,y,'maxima');
%%% Step2: Finding associated minima or maxima of each contour
bnd = struct('xc',[],'yc',[],'cval',[],'xp',[],'yp',[],'valp',[]);
for ii=1:numel(ct_x)
    ind_in = find( inpolygon(xp,yp,ct_x{ii},ct_y{ii}) );
    [~,indmax] = max(valp(ind_in));
    indmax = ind_in(indmax);
    if isempty(indmax); continue; end;
    bnd.xc{end+1} = ct_x{ii};
    bnd.yc{end+1} = ct_y{ii};
    bnd.cval(end+1,1) = ct_val{ii};
    bnd.xp(end+1,1) = xp(indmax);
    bnd.yp(end+1,1) = yp(indmax);
    bnd.valp(end+1,1) = valp(indmax);
end

