%
%  File    :   pca_test.m 
%  Author  :   Martin R. Cole
%  Date    :   19/04/2024
%  

% Add enCoRe to the Path
addpath('interpolation/')
addpath('kde/')
addpath('simulations/')
addpath('tangent_basis')
addpath('utils/')

%% Global settings

input_folder = '../grid';
input_folder2 = '../subject_ids';
input_folder3 = '../example';

% CHANGE THIS TO YOUR OWN WORKING DIRECTORY
output_folder = '../results';

% grid to do evaluations over
[lh_grid,~] = read_vtk(sprintf('%s/lh_grid_avg_0.94.vtk', input_folder));
[rh_grid,~] = read_vtk(sprintf('%s/rh_grid_avg_0.94.vtk', input_folder));

% load all the subject IDs
sublist = load(sprintf('%s/hcp_all_ids.csv', input_folder2));

% using the same subjects the we did for the
% template and registration in the previous scripts
N_subs = 10;

%rng(54569783); % random seed for reproducibility
%selected_subs = sublist(randsample(length(sublist), N_subs));
selected_subs = sublist(1:N_subs); % the example data

%% Setup (same as for the template)

% normalise the grid to radius 1
lh_mesh.V = normr(lh_grid.vtx');
lh_mesh.T = lh_grid.tri;
rh_mesh.V = normr(rh_grid.vtx');
rh_mesh.T = rh_grid.tri;

% initialise the Concon object
concon = Concon(SphericalGrid(lh_mesh,30), ...
                SphericalGrid(rh_mesh,30), 1e-10);

P = length(lh_mesh.V) + length(rh_mesh.V);
N_elems = ((P-1)*P) / 2;

elem_idx = logical(triu(ones(P,P), 1));

%% PCA for original SC

clear PCA_mat
PCA_mat = zeros(N_subs,N_elems);

for i = 1:N_subs
    tmp = load(sprintf('%s/SC_%s.mat', input_folder3, string(selected_subs(i))));

    % load SC the same way as with the template
    SC = (tmp.SC + tmp.SC') - 2*diag(diag(tmp.SC));
    SC = SC ./ sum(SC, 'all');
    %SC = log(SC + 1);

    PCA_mat(i,:) = SC(elem_idx);
    disp(i)
end

[~,orig_score,~,orig_explained,~] = fastpca(PCA_mat);

%% PCA for registered SC

clear PCA_mat warp_mat
PCA_mat = zeros(N_subs,N_elems);
warp_mat = zeros(N_subs,P);

for i = 1:N_subs
    tmp = load(sprintf('%s/SC_%s.mat', input_folder3, string(selected_subs(i))));
    load(sprintf('%s/registered_warp_%s.mat', output_folder, string(selected_subs(i))),'lh_warp','rh_warp');

    % load SC the same way as with the template
    SC = (tmp.SC + tmp.SC') - 2*diag(diag(tmp.SC));
    SC = SC ./ sum(SC, 'all');
    %SC = log(SC + 1);     
    
    % register the SC with the saved warp
    SC = concon.evaluate(SC,lh_warp,rh_warp);

    % calculate geodesic (great circle) distance between
    % the original mesh and the warped mesh vertices
    warp_mat(i,1:length(lh_mesh.V)) = acos(max(-1, min(1, dot(lh_warp.V,lh_mesh.V,2))));
    warp_mat(i,(length(lh_mesh.V)+1):end) = acos(max(-1, min(1, dot(rh_warp.V,rh_mesh.V,2))));
  
    PCA_mat(i,:) = SC(elem_idx);
    disp(i)
end

[~,reg_score,~,reg_explained,~] = fastpca(PCA_mat);
[~,warp_score,~,~,warp_explained] = pca(warp_mat);
