%% Example script to generate, illustrate, and compute mutual admittances


%% Physical constants
f = 10E9;
mu = 1.25663706212*1E-6;
epsilon  = 8.8541878128*1E-12;
lambda = 1/sqrt(mu*epsilon) / f;
k = 2*pi*f*sqrt(epsilon*mu);
omega = 2*pi*f;
c = 1/sqrt(mu*epsilon);
eta = sqrt(mu/epsilon);
Y_intrinsic_source = 35.3387;

l_m = 1; % Magnetic length of dipoles. This is just a normalization constant.

%% Define system

a = 0.7318 * lambda;
b = 0.1668 * lambda;

nRf = 2; % number of waveguides
sRf = lambda;  % spaceing of waveguides
nDMA = 5; % per waveguide (does not have to be the same for all waveguides in general).
sDMA = 0.5*lambda; % DMA element spacing along waveguide
l = sDMA * (nDMA+2);


x_dma = (1:nDMA) * sDMA; x_dma = x_dma - mean(x_dma) + l/2;
y_dma = b;
z_dma = (0:(nRf-1)) * sRf + a/2;

[Xt,Yt,Zt] = ndgrid(x_dma, y_dma, z_dma);
xyz_dma = [Xt(:), Yt(:), Zt(:)];

x_tx = 0;
y_tx = b/2;
z_tx = (0:(nRf-1)) * sRf + a/2;

[Xt,Yt,Zt] = ndgrid(x_tx, y_tx, z_tx);
xyz_tx = [Xt(:), Yt(:), Zt(:)];


%% Illustrate DMA topology

figure; hold on; legend;
xlabel('x'); ylabel('y'); zlabel('z');
scatter3(xyz_dma(:,1),xyz_dma(:,2),xyz_dma(:,3), 'displayName','DMA elements')

scatter3(xyz_tx(:,1),xyz_tx(:,2),xyz_tx(:,3), 'displayName','Tx')

% Draw the waveguides
xtot = zeros(0,4);
ytot = zeros(0,4);
ztot = zeros(0,4);
x_max = l;
for idn = 1:nRf
    zt = [0, a, a, 0; ...
        0, a, a, 0; ...
        0, 0, 0, 0; ...
        a, a, a ,a] + sRf*(idn-1);
    yt = [0, 0, 0, 0; ...
        b, b, b, b; ...
        0, b, b, 0; ...
        0, b, b, 0];
    xt = [0, 0, x_max, x_max; ...
        0, 0, x_max, x_max;...
        0, 0, x_max, x_max;...
        0, 0, x_max, x_max];
    xtot = [xtot;xt];
    ytot = [ytot;yt];
    ztot = [ztot;zt];
end
h1 = patch('XData',xtot','YData',ytot','ZData',ztot','displayName','waveguides');
h1.FaceAlpha = 0.2;
view([-21,37]);

axis equal

%% Compute mutual admittances

[Y_tt, Y_st, Y_ss] = DMA_admittance(f, a, b, l_m, l, xyz_dma, xyz_tx);


%% Generate rayleigh channel to two users
% Note that sigma_alpha^2 is chosen as to normalize the variance.

nRx = 2;
nTx = nRf * nDMA;
H_uncorr = (randn(nRx, nTx) + 1i * randn(nRx, nTx)) .* sqrt(1/2);

if nTx == 1
    sq_Sigma = 1;
else
    dZ = squareform(pdist(xyz_tx(:,3)));
    dR = squareform(pdist(xyz_tx));
    k = 2*pi/lambda;
    Sigma = 3./2.*((1 + (-k.^2.*dZ.^2 - 1)./(dR.^2.*k.^2) + 3.*dZ.^2./(dR.^4.*k.^2)).*sin(k.*dR)./(k.*dR) ...
        + cos(k.*dR).*(1./(k.*dR) - 3.*dZ.^2./(dR.^3.*k))./(k.*dR)); % Normalized
    
    Sigma(isnan(Sigma)) = 1;
    sq_Sigma = real(Sigma^(1/2)); % Real operator due to imaginary part being a product of quantization
end
H = sq_Sigma * H_uncorr;

Y_rs = H;
Y_rr = eye(nRx) * k * 2*pi*f * epsilon / 6*pi;
Y_r = eye(nRx) .* Y_rr'; % Complex conjugate match

% Equivalent channel from the RF chain to the users:
Y_s = 0 + 1i*diag(randn(nTx,1) * 10); % Random configuration
H_eq = (Y_r + Y_rr) \ ( (Y_rs / (Y_s + Y_ss)) * Y_st); 

% Computing approximate reflection coefficient assuming no cross-waveguide
% coupling
Y_p = Y_tt - (Y_st.' / (Y_s + Y_ss)) * Y_st;
Y_in = eye(nRf) .* Y_p;
Gamma = (Y_in - eye(nRf)*Y_intrinsic_source) / (Y_in + eye(nRf)*Y_intrinsic_source);

% Choosing random beam
B = randn(nRf,nRx) + 1i*randn(nRf,nRx);

% Computing received, transmit, and supplied power
x = [1;1];
y = H_eq * B * x;

P_r = 1/2 * real(Y_r) * abs(y).^2
P_t = 1/2 * real(x' * B' * Y_p * B * x)
P_s = 1/2 * real(x' * B' * ((eye(nRf) - Gamma' * Gamma) \ Y_p) * B * x)
