function [Y_tt, Y_st, Y_ss] =... 
    DMA_admittance(f, a, b, l, S_mu, xyz_dma, xyz_rf)
% Computes the mutual admittances in the DMA topology where all actors are
% modelled as magnetic dipoles
%
% RJW 03/02/2022 - Added length variable for normalization purposes.
%

%% Physical constants - Dont edit
mu = 1.25663706212*1E-6;
epsilon  = 8.8541878128*1E-12;
k = 2*pi*f*sqrt(epsilon*mu);

kx = sqrt(k^2 - (pi/a)^2);    % only for TE_10

%% Y_tt calculation
% Position of RF chain in the width of the waveguide (0 <= z_rf <= a)
z_rf = a/2;

y_tt = -l^2 * 2i*kx*sin(pi/a*z_rf)^2*cos(kx*S_mu) / (a*b*2*pi*f*mu*sin(kx*S_mu));

Y_tt = diag(y_tt * ones(size(xyz_rf,1),1));

%% Y_st calculation 
Y_st = zeros(size(xyz_dma,1),size(xyz_rf,1));

% r, rhat are vectors such that 0 <= x, xhat <= S_mu
%                               0 <= y, yhat <= b
%                               0 <= z, zhat <= a
Gw_e2 = @(r,rhat) -kx*sin(pi/a*rhat(3))*sin(pi/a*r(3))*...
                      (cos(kx*(rhat(1) + r(1)-S_mu)) + ...
                       cos(kx*(S_mu-abs(rhat(1) - r(1))))) / (a*b*k^2*sin(kx*S_mu)); 

% We compute the normalized x coordinate of the antennas wrt their RF chain
N_ant_wg = size(xyz_dma,1)/size(xyz_rf,1);
temp = zeros(size(xyz_dma,1),1);
temp(1:N_ant_wg:end) = xyz_rf(:,1);
x_ant_norm = xyz_dma(:,1) - ...
    filter(ones(N_ant_wg,1), 1, temp);
                   
% We assume all the waveguides are equal, i.e., same position for the
% elements and RF chain, so that Y_st is a block matrix with same elements
% within each block
for row = 1:N_ant_wg
    y_val = l^2 * 1i*2*pi*f*epsilon*Gw_e2([x_ant_norm(row), 0, z_rf], [0, 0, z_rf]);
    for k_rf = 1:size(xyz_rf,1)
        Y_st(row + (k_rf-1)*N_ant_wg,k_rf) = y_val;
    end
end

%% Y_ss calculation 
Ga_e2 = @(r,rhat) ((norm(r-rhat)^2 - (r(3)-rhat(3))^2)/norm(r-rhat)^2 - ...
                    1i*(norm(r-rhat)^2 - 3*(r(3)-rhat(3))^2)/(norm(r-rhat)^3*k)...
                    - (norm(r-rhat)^2 - 3*(r(3)-rhat(3))^2)/(norm(r-rhat)^4*k^2))...
                    *exp(-1i*k*norm(r-rhat))/(4*pi*norm(r-rhat));
                
Y_ss = zeros(size(xyz_dma,1),size(xyz_dma,1));

for row = 1:size(xyz_dma,1)
   for col = 1:size(xyz_dma,1)
      % If diagonal element -> self-admittance
      if(row == col)
          Y_ss(row,col) = l^2 * k*2*pi*f*epsilon/(3*pi) + ...
      1i*2*pi*f*epsilon*Gw_e2([xyz_dma(row,1), 0, z_rf], [xyz_dma(col,1), 0, z_rf]);
          
      % If different waveguide
      elseif (abs(xyz_dma(row,1) - xyz_dma(col,1))>= S_mu || ...
              abs(xyz_dma(row,3) - xyz_dma(col,3))>= b/2)
          
          Y_ss(row,col) = l^2 * 1i*2*pi*f*epsilon*...
              2*Ga_e2(xyz_dma(row,:),xyz_dma(col,:));
          
      else
          Y_ss(row,col) = l^2 * 1i*2*pi*f*epsilon*...
            (2*Ga_e2(xyz_dma(row,:),xyz_dma(col,:)) + ...
            Gw_e2([xyz_dma(row,1), 0, z_rf], [xyz_dma(col,1), 0, z_rf]));
      end
   end
end
                
end