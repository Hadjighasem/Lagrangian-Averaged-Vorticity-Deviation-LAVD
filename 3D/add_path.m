%------------------------------------Set path
function add_path

fp = mfilename('fullpath');
rootdir = fileparts(fp);
p{1} = fullfile(rootdir,'data');
p{2} = fullfile(rootdir,'doc');
p{3} = fullfile(rootdir,'examples');
p{4} = fullfile(rootdir,'src');

for i = 1:4
    addpath(rootdir,p{i});
end
%------------------------------------
fprintf('----------------------------*---------------------------------\n');
fprintf('data,...doc,...examples,...src,...All paths have been added');
fprintf('\n')
fprintf('----------------------------*---------------------------------\n');

end