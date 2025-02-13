%
%  File    :   pop_template.m 
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

% select N random subjects to create the template
N_subs = 10;

%rng(54569783); % random seed for reproducibility
%selected_subs = sublist(randsample(length(sublist), N_subs));
selected_subs = sublist(1:N_subs); % the example data

%% Setup

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

% preallocate array to fill (more memory efficient)
P = length(lh_mesh.V) + length(rh_mesh.V);
SC = zeros(P,P,N_subs);

%% Load Data

for i = 1:N_subs
    tmp = load(sprintf('%s/SC_%s.mat', input_folder3, string(selected_subs(i))));

    % make the matrix symmetric and remove diagonal
    SC(:,:,i) = (tmp.SC + tmp.SC') - 2*diag(diag(tmp.SC));
    SC(:,:,i) = SC(:,:,i) ./ sum(SC(:,:,i), 'all');
    
    % perform log transformation to make SC smoother (may not be needed)
    %SC(:,:,i) = log(SC(:,:,i) + 1);

    disp(i);
end

%% Generate Template

% create a template using a maximum of 200 iterations
template = encore.get_template(SC,200);
save(sprintf('%s/template.mat', output_folder),'template','-v7.3');
