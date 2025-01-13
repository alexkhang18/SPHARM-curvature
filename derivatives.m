function [r_phi,r_phiphi,r_theta,r_thetatheta,r_phitheta] = derivatives(Ymn, amn, SPHARM_verts,deg)

% Convert Cartesian parameterized coordinates to spherical coordinates
[phi, theta] = cart2sph(SPHARM_verts(:,1), SPHARM_verts(:,2), SPHARM_verts(:,3));
ind = find(phi<0);
phi(ind) = phi(ind) + 2*pi;
theta = pi/2 - theta;

% computes SPHARM basis functions
function v = get_Y(l,m,n)
    p = get_P(l,m,n);
    v = get_f(l,m,n) * p * exp(1i*m*phi(n));
end

% computes Legendre polynomial portion of SPHARM basis functions
function v = get_P(l,m,n)
    Plm = legendre(l,cos(theta(n)));
    if m >= 0
        v = Plm(m+1);
    end
    if m<0 && mod(m,2) == 1
        v = Plm(abs(m)+1) * (-1) * factorial(l-abs(m)) / factorial(l+abs(m));  
    end
    
    if m<0 && mod(m,2) == 0
        v = Plm(abs(m)+1) * factorial(l-abs(m)) / factorial(l+abs(m)); 
    end
end

% computes non-Legendre polynomial portion of SPHARM basis functions
function v = get_f(l,m,~)
    a = (2*l+1)*factorial(l-m);
    b = 4*pi*factorial(l+m);
    v = sqrt(a/b);
end

% initiate derivatives
r = zeros(length(theta),3);
r_phi = zeros(length(theta),3);
r_theta = zeros(length(theta),3);
r_phiphi = zeros(length(theta),3);
r_thetatheta = zeros(length(theta),3);
r_phitheta = zeros(length(theta),3);

% initiate temporary derivative values
N = (deg+1)^2;
tmp_r_1 = zeros(1,N);tmp_r_2 = zeros(1,N);tmp_r_3 = zeros(1,N);
tmp_r_phi_1 = zeros(1,N);tmp_r_phi_2 = zeros(1,N);tmp_r_phi_3 = zeros(1,N);
tmp_r_phiphi_1 = zeros(1,N);tmp_r_phiphi_2 = zeros(1,N);tmp_r_phiphi_3 = zeros(1,N);
tmp_r_theta_1 = zeros(1,N);tmp_r_theta_2 = zeros(1,N);tmp_r_theta_3 = zeros(1,N);
tmp_r_thetatheta_1 = zeros(1,N);tmp_r_thetatheta_2 = zeros(1,N);tmp_r_thetatheta_3 = zeros(1,N);
tmp_r_phitheta_1 = zeros(1,N);tmp_r_phitheta_2 = zeros(1,N);tmp_r_phitheta_3 = zeros(1,N);

% loop to compute derivatives for every nodal point (theta), every degree (l), and
% every order (m)
for n=1:length(theta) % iterates for every point

    ii = 1;

    for l=0:deg % iterates for every degree

        for m=-l:l % iterates for every order
            
            % concatenates coefficients
            amn_1 = amn(ii,1);
            amn_2 = amn(ii,2);
            amn_3 = amn(ii,3);
            ymn = Ymn(n,ii);

            % Computes components of derivatives
            A = (l+1)*cos(theta(n))*get_P(l,m,n)-(l-m+1)*get_P(l+1,m,n);
            B = (l+1+(l+1)^2*cos(theta(n))^2)*get_P(l,m,n)-2*cos(theta(n))*...
                (l-m+1)*(l+2)*get_P(l+1,m,n)+(l-m+1)*(l-m+2)*get_P(l+2,m,n);
            E = exp(1i*m*phi(n));
            F = get_f(l,m,n);
            S = sin(theta(n));       

            % Computes derivatives for specific node, degree, and order and saves them
            tmp_r_1(ii) = real(amn_1 * ymn);
            tmp_r_2(ii) = real(amn_2 * ymn);
            tmp_r_3(ii) = real(amn_3 * ymn);
            tmp_r_phi_1(ii) = 1i * m * amn_1 * ymn;
            tmp_r_phi_2(ii) = 1i * m * amn_2 * ymn;
            tmp_r_phi_3(ii) = 1i * m * amn_3 * ymn;
            tmp_r_phiphi_1(ii) = - m^2 * amn_1 * ymn;
            tmp_r_phiphi_2(ii) = - m^2 * amn_2 * ymn;
            tmp_r_phiphi_3(ii) = - m^2 * amn_3 * ymn;
            tmp_r_theta_1(ii) = - F / S * amn_1 * A * E;
            tmp_r_theta_2(ii) = - F / S * amn_2 * A * E;
            tmp_r_theta_3(ii) = - F / S * amn_3 * A * E;
            tmp_r_thetatheta_1(ii) = F / S^2 * amn_1 * B * E;
            tmp_r_thetatheta_2(ii) = F / S^2 * amn_2 * B * E;
            tmp_r_thetatheta_3(ii) = F / S^2 * amn_3 * B * E;
            tmp_r_phitheta_1(ii) = - 1i * m * F / S * amn_1 * A * E;
            tmp_r_phitheta_2(ii) = - 1i * m * F / S * amn_2 * A * E;
            tmp_r_phitheta_3(ii) = - 1i * m * F / S * amn_3 * A * E;

            ii = ii+1;
        end

    end

    % sums up all derivatives
    r (n,1) = sum(tmp_r_1);r (n,2) = sum(tmp_r_2);r (n,3) = sum(tmp_r_3);
    r_phi(n,1) = sum(tmp_r_phi_1);r_phi(n,2) = sum(tmp_r_phi_2);r_phi(n,3) = sum(tmp_r_phi_3);
    r_phiphi(n,1) = sum(tmp_r_phiphi_1);r_phiphi(n,2) = sum(tmp_r_phiphi_2);r_phiphi(n,3) = sum(tmp_r_phiphi_3);
    r_theta(n,1) = sum(tmp_r_theta_1);r_theta(n,2) = sum(tmp_r_theta_2);r_theta(n,3) = sum(tmp_r_theta_3);
    r_thetatheta(n,1) = sum(tmp_r_thetatheta_1);r_thetatheta(n,2) = sum(tmp_r_thetatheta_2);r_thetatheta(n,3) = sum(tmp_r_thetatheta_3);
    r_phitheta(n,1) = sum(tmp_r_phitheta_1);r_phitheta(n,2) = sum(tmp_r_phitheta_2);r_phitheta(n,3) = sum(tmp_r_phitheta_3);

end

end