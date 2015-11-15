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
bnd = struct('xc',[],'yc',[],'cval',[],'xp',[],'yp',[],'valp',[]);
for k=1:numel(s)
    % check if the contour is closed
    if ( s(k).x(1)==s(k).x(end) ) && ( s(k).y(1)==s(k).y(end) ) 
        Length = sum( sqrt(diff(s(k).x,1).^2+diff(s(k).y,1).^2) );
        if Length>MinLength
            % Check if the contour is convex
            if IsContourConvex(s(k).x,s(k).y,DeficiencyThresh) 
                bnd.xc{end+1} = s(k).x';
                bnd.yc{end+1} = s(k).y'; 
                bnd.cval(end+1,1) = s(k).v';
            end
        end
    end
end
Nct_filt1 = numel(bnd.xc);
disp(sprintf('... %3d contours are closed and convex.',Nct_filt1));
%% Step2: Finding local maxima of the LAVD field.
[xp,yp,valp] = Find2DPeak(VMatrix,x,y,'maxima');
%- Keeping closed contours that encompass only "1" local maximum.
InMax = cellfun(@(x,y) inpolygon(xp,yp,x,y),bnd.xc,bnd.yc,'UniformoutPut',false);
InMax = cat(1,InMax{:});       % rows --> closed contours & columns --> local maxima
N_InMax = sum(InMax,2);        % Number of local maxima in each contour
indEliminate_1 = N_InMax~=1;
bnd.xc(indEliminate_1)   = [];
bnd.yc(indEliminate_1)   = [];
bnd.cval(indEliminate_1) = [];

[indmax,~] = find( InMax(~indEliminate_1,:)'==1 );
bnd.xp = xp(indmax);
bnd.yp = yp(indmax);
bnd.valp = valp(indmax);
%% Step3: Selecting the outermost contour for each nested family of contours
%- selecting the first point of each contour as a query point
Nct_filt2 = numel(bnd.xc);
if Nct_filt2~=0
    xq = cellfun(@(x) x(1),bnd.xc,'UniformoutPut',true);  
    yq = cellfun(@(x) x(1),bnd.yc,'UniformoutPut',true);  
end
%%
indEliminate_2 = false(Nct_filt2,1);
for ii=1:Nct_filt2
    in = inpolygon(xq,yq,bnd.xc{ii},bnd.yc{ii});  in(ii)=0;
    indEliminate_2(in) = 1;
end
bnd.xc(indEliminate_2)   = [];
bnd.yc(indEliminate_2)   = [];
bnd.cval(indEliminate_2) = [];

bnd.xp(indEliminate_2)   = [];
bnd.yp(indEliminate_2)   = [];
bnd.valp(indEliminate_2) = [];

end
