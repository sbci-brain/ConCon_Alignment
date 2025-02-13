% MATLAB R2022b
%
% Interpolation SBCI to Full mesh resolution - Performs barycentric
% interpolation to upsample data from SBCI resolution to FreeSurfer
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 
% CITE REFERENCES HERE
% 
% Author: M.R. Cole, PhD
% Department of Psychology, University of Rochester
% email: martin_cole@urmc.rochester.edu
% 30-Apr-2024 ; Last revision: 02-May-2024 

% Load low resolution SBCI mesh grid
input_folder = '../grid';

[lh_grid,~] = read_vtk(sprintf('%s/lh_grid_avg_0.94.vtk', input_folder));
[rh_grid,~] = read_vtk(sprintf('%s/rh_grid_avg_0.94.vtk', input_folder));

% Convert struct to format used by AABB tree
lh_grid.V = lh_grid.vtx';
rh_grid.V = rh_grid.vtx';
lh_grid.T = lh_grid.tri;
rh_grid.T = rh_grid.tri;

% Create AABB trees for interpolation
lh_tree = AABBtree(lh_grid);
rh_tree = AABBtree(rh_grid);

% load required SBCI data for mapping and analysis
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('/overflow/zzhanglab/encore/SBCI_Toolkit/example_data', '0.94');
sbci_surface = load_sbci_surface('../SBCI_Toolkit/example_data');

lh_points = normr(sbci_surface.sphere.lh_surf.vtx.');
rh_points = normr(sbci_surface.sphere.rh_surf.vtx.');

% Perform upsampling
example_data = load('../SBCI_Toolkit/example_data/example_sfc.mat', 'sfc');

% create a mask for the corpus callosum
mask = (sbci_parc(8).labels == 1) | (sbci_parc(8).labels == 9);
example_data.sfc(mask) = NaN;

[Vq,Tq] = lh_tree.get_barycentric_data(lh_points, false);

lh_new_data = Vq(:,1) .* example_data.sfc(Tq(:,1)) + ...
              Vq(:,2) .* example_data.sfc(Tq(:,2)) + ...
              Vq(:,3) .* example_data.sfc(Tq(:,3));

[Vq,Tq] = rh_tree.get_barycentric_data(rh_points, false);
rh_new_data = Vq(:,1) .* example_data.sfc(Tq(:,1)+size(lh_grid.V,1)) + ...
              Vq(:,2) .* example_data.sfc(Tq(:,2)+size(lh_grid.V,1)) + ...
              Vq(:,3) .* example_data.sfc(Tq(:,3)+size(lh_grid.V,1));

% plot upsampled data
plot_cortical_highres(sbci_surface.inflated,lh_new_data,rh_new_data, 'figid', 1, 'cmap', 'parula')

% plot original data using old upsampling
plot_cortical(sbci_surface.inflated, sbci_mapping, example_data.sfc, 'figid', 2, 'cmap', 'parula')