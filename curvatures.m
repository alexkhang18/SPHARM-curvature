function [mean_k,gauss_k] = curvatures(r_theta,r_phi,r_thetatheta,r_phiphi,r_phitheta)

% take real portion
r_phi = real(r_phi);
r_theta = real(r_theta);
r_phiphi = real(r_phiphi);
r_thetatheta = real(r_thetatheta);
r_phitheta = real(r_phitheta);

% calculate derivative of normal with respect to theta and phi
xi = cross(r_theta,r_phi,2);
xi_theta = cross(r_thetatheta,r_phi,2) + cross(r_theta,r_phitheta);
xi_phi = cross(r_phitheta,r_phi,2) + cross(r_theta,r_phiphi,2);

nn = xi ./ sqrt(dot(xi,xi,2));

top1 = xi_theta .* sqrt(dot(xi,xi,2));
top2 = dot(xi_theta,xi,2) .* nn;
bottom1 = dot(xi,xi,2);
nn_theta = (top1 - top2)./bottom1;

top1 = xi_phi .* sqrt(dot(xi,xi,2));
top2 = dot(xi_phi,xi,2) .* nn;
bottom1 = dot(xi,xi,2);
nn_phi = (top1 - top2) ./ bottom1;

% flip normal directions
nn = -nn;
nn_theta = -nn_theta;
nn_phi = -nn_phi;

% First fundamental form
E = dot(r_theta,r_theta,2);
F = dot(r_theta,r_phi,2);
G = dot(r_phi,r_phi,2);

% Second fundamental form
e = - dot(r_theta,nn_theta,2);
f = 0.5 .* (dot(r_theta,nn_phi,2) + dot(r_phi,nn_theta,2));
g = - dot(r_phi,nn_phi,2);

% Computes components of shape operator
S11 = (E.*G-F.^2).^(-1).*[e.*G-f.*F];
S12 = (E.*G-F.^2).^(-1).*[f.*G-g.*F];
S21 = (E.*G-F.^2).^(-1).*[f.*E-e.*F];
S22 = (E.*G-F.^2).^(-1).*[g.*E-f.*F];

% Computes curvature metrics
for j = 1:length(E)
    S_temp = [S11(j),S12(j);S21(j),S22(j)]; % assembles shape operator
    mean_k(j) = 1/2*trace(S_temp); % computes mean k
    gauss_k(j) = det(S_temp); % computes gauss k
end

end