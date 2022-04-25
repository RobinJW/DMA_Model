%% Parameters
f = 10e9;                        % Frequency of operation
c = 299792458;                   % Light speed in vacuum
mu = 4*pi*1E-7;                  % Vacuum permeability
epsilon = 1/ (c^2 * mu);         % Vacuum permittivity
lambda = c/f;                    % Wavelength
k = 2*pi/lambda;                 % Wavenumber
a = 0.73*lambda;                 % Width of waveguides (only TE_10 mode)
b = 0.17*lambda;                 % Height of waveguides (only TE_10 mode)
channel_type = 1;                % Type of channel:
                                 % 0 -> LoS 
                                 % 1 -> Rayleigh
N = 4;                           % Number of RF chains / waveguides
Lmu = 8;                         % Number of elements per waveguide
wvg_spacing = lambda;            % Spacing between waveguides
elem_spacing = lambda/2;         % Spacing between the elements
l = 1;                           % Length of dipoles -> just normalization
M = 3;                           % Number of static users
Plot_topology = 1;               % Boolean to plot the chosen setup

Y_s = diag(1i*randn(N*Lmu,1));   % Load admittances of DMA element
                                 % Has to be a diagonal matrix of 
                                 % N*Lmu x N*Lmu (total number of elements)

Y_intrinsic_source = 35.3387;    % Intrinsic impedance of source 
                                 % matched to waveguide of width a  =
                                 % 0.73*lambda and height b = 0.17*lambda


%% DMA and users coordinates
site_xyz = [0 0 10];             % [x y z] coordinates of bottom right 
                                 % corner of DMA
S_mu = (Lmu+1)*elem_spacing;     % Length of waveguides

% Coordinates of DMA elements and RF chains
[ant_xyz, rf_xyz] = Topologies_DMA(site_xyz,N, Lmu, wvg_spacing,...
                elem_spacing, S_mu, a, b, Plot_topology);
            
% Users positions (In this example, they are set randomly)
x_lim = [-20 20];
y_lim = [20 60];
user_xyz = [x_lim(1)+(x_lim(2)-x_lim(1))*rand(M,1) ...
            y_lim(1)+(y_lim(2)-y_lim(1))*rand(M,1) 1.5*ones(M,1)]; 

%% Calculation of Admittances 

% Calculating Y_tt, Y_st and Y_ss according to Eqs. (35)-(42)
[Y_tt, Y_st, Y_ss] = DMA_admittance(f, a, b, l, S_mu, ant_xyz, ...
                                    rf_xyz, mu, epsilon);
                                
% Calculating Y_rr according to Eqs. (44)-(46)
Y_rr = Coupling_Dipoles(f, l, user_xyz, mu, epsilon);
                                
% Choosing Y_r (load admittance of users) as conjugate of self-admittance
Y_r = Y_rr'.*eye(M);

% Calculation of Y_rs (Wireless channel)
Y_rs = GenChannel(channel_type, lambda, ant_xyz, user_xyz);

%% Equivalent channel according to Eq. (60)
Heq = eye(M)/(Y_r + Y_rr) * (Y_rs/(Y_s + Y_ss)*Y_st);

%% Computing received, transmitted and supplied power:

% Computing approximate reflection coefficient assuming no cross-waveguide
% coupling
Y_p = Y_tt - (Y_st.' / (Y_s + Y_ss)) * Y_st;
Y_in = eye(N) .* Y_p;
Gamma = (Y_in - eye(N)*Y_intrinsic_source) / (Y_in + eye(N)*Y_intrinsic_source);

% Choosing random beam
B = randn(N,M) + 1i*randn(N,M);

% Computing received, transmit, and supplied power
x = randn(M,1);
y = Heq * B * x;

P_r = 1/2 * real(Y_r) * abs(y).^2
P_t = 1/2 * real(x' * B' * Y_p * B * x)
P_s = 1/2 * real(x' * B' * ((eye(N) - Gamma' * Gamma) \ Y_p) * B * x)
