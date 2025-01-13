%% SPHARM + Curvature

clc
clear all
close all

%% set SPHARM degree
deg = 10; % a higher degree will result in a better fit

%% reads in surface mesh

% read mesh
mesh = readSurfaceMesh('T_1.ply');

% vertices or nodes
vertices = mesh.Vertices;

% faces or connectivity
faces = double(mesh.Faces);

% centers the mesh
vertices_mean = mean(vertices);
vertices = vertices - vertices_mean;

% writes centered mesh to a .vtk file (open in Paraview)
writeVTK('original_shape.vtk', vertices, faces)

%% Spherical paramterization via Allen Institute Method
% Reference - Viana, M. P., et al., Nature 613, 345–354 (2023).

% converts mesh vertices from cartesian to spherical coordinates 
[PHI,THETA,R] = cart2sph(vertices(:,1),vertices(:,2),vertices(:,3));

% generates sphere points for interpolation
PHIq = -pi:pi/64:pi;   
THETAq = -pi/2:pi/64:pi/2;   
PHIq(end) = []; % removes last point
THETAq(end) = []; % removes last point
[Phi_new, Theta_new] = meshgrid(PHIq, THETAq); % creates grid
theta_grid = Theta_new(:); % turns grid into list of points
phi_grid = Phi_new(:); % turns grid into list of points
phi_grid = [phi_grid;0;0]; % adds in the north and south poles
theta_grid = [theta_grid;-pi/2;pi/2]; % adds in the north and south poles

% interpolatation onto sphere
F_interp = scatteredInterpolant(THETA, PHI, R, 'nearest', 'nearest');
R_grid = F_interp(theta_grid, phi_grid);

% converts back to catesian points
[x,y,z] = sph2cart(phi_grid,theta_grid,ones(size(phi_grid)));
sphere_verts = [x,y,z];

% slightly shifts the poles so that normals can be computed without error
indx1 = find(sphere_verts(:,3)==1);
indx_1 = find(sphere_verts(:,3)==-1);
sphere_verts(indx1,:) = normr(sphere_verts(indx1,:) - 1e-6);
sphere_verts(indx_1,:) = normr(sphere_verts(indx_1,:) - 1e-6);
[PHI,THETA,R] = cart2sph(sphere_verts(:,1),sphere_verts(:,2),sphere_verts(:,3));

% creates faces or connectivity for new sphere
fMat = boundary(sphere_verts,0);

% creates new interpolated shape
[x,y,z] = sph2cart(PHI,THETA,R_grid);
sphere_verts = [x,y,z];    

%% Compute SPHARM basis functions and coefficients
% Reference - Shen, L., et al., Evolution 63, 1003–1016 (2009).

Zq = calculate_SPHARM_basis(sphere_verts, deg); % function from Shen, L. et al.

coeffs = Zq\sphere_verts;

SPHARM_verts = real(Zq*coeffs);

%% Computes curvatures 
% Reference - Garboczi, E. J., Cem. Concr. Res. 32, 1621–1638 (2002).
  
% Calculate derivatives
[r_phi,r_phiphi,r_theta,r_thetatheta,r_phitheta] = derivatives(Zq, coeffs, sphere_verts, deg);

% Calculate local mean curvature H and Gaussian curvature K
[mean_k,gauss_k] = curvatures(r_theta,r_phi,r_thetatheta,r_phiphi,r_phitheta);

%% Outputs to .vtk to be viewed in Paraview

writeVTK('SPHARM_curvatures_degree_10.vtk', SPHARM_verts, fMat,...
'mean_k',mean_k',...
'gauss_k',gauss_k')

   

    

