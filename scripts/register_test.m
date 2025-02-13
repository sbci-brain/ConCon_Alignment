%
%  File    :   register_test.m 
%  Author  :   Martin R. Cole
%  Date    :   19/04/2024
%  

% Add enCoRe to the Path
addpath('interpolation/')
addpath('kde/')
addpath('simulations/')
addpath('tangent_basis')
addpath('utils/')
addpath('./')

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

% Can select any subjects to register to the pregenerated template,
% but for this test, we'll use the same subjects the we did for the
% template in the previous script
N_subs = 10;

rng(54569783); % random seed for reproducibility
selected_subs = sublist(randsample(length(sublist), N_subs));

%% Setup (same as for the template)

% normalise the grid to radius 1
lh_mesh.V = normr(lh_grid.vtx');
lh_mesh.T = lh_grid.tri;
rh_mesh.V = normr(rh_grid.vtx');
rh_mesh.T = rh_grid.tri;

% initialise the enCoRe object:
% spherical harmonics to order = 30
% step size for gradient descent = 0.01
% maximum iterations = 1000
% tolerance for convergence = 0.000006
encore = Encore(lh_mesh,rh_mesh,30,0.005,1000,1e-6);

%% Register Population to the Template

load(sprintf('%s/template.mat', output_folder),'template');

for i = 1:N_subs
    tmp = load(sprintf('%s/SC_%s.mat', input_folder3, string(selected_subs(i))));

    % load SC the same way as with the template
    SC = (tmp.SC + tmp.SC') - 2*diag(diag(tmp.SC));
    SC = SC ./ sum(SC, 'all');
    %SC = log(SC + 1);

    disp(i)
    tic    
    [~,lh_warp,rh_warp,~] = encore.register(template,SC,'is_template',true);
    toc

    save(sprintf('%s/registered_warp_%s.mat', output_folder, string(selected_subs(i))), 'lh_warp', 'rh_warp', '-v7.3');    
end
