%% Script to generate figures for the paper.
%
% NOTE here that the coordinate system is rotated 90 degrees in the 
% xz-plane to align with the CST simulation used for comparison.
%
% To clarify, the DMA is positioned in the xz-plane such that:
% Waveguide width is along the x-coordinate
% Waveguide height is along the y-coordinate
% Waveguide length is along the z-coordinate

plot_stuff = true;
remove_spherical_guides = 1;
remove_cartesian_guides = 0;

addpath('simulationData');
setFigOptions;
colours = [   0         0    1.0000 ;
    1.0000         0         0;
    0    1.0000         0;
    0         0    0.1724];
colourIndx = 0;
lineSpec = {'--','-',':','-.'};
markerSpec = {'+','','o','*','x','s','d','^','v','>','<','p','h'};
right = 1; left = 0;

linewidth = 262*0.03514;

Y_intrinsic_source = 35.3387;

orderM = 1;
orderN = 0;


%% Physical constants
f = 10E9;
mu = 1.25663706212*1E-6;
eps  = 8.8541878128*1E-12;
lambda = 1/sqrt(mu*eps) / f;
k = 2*pi*f*sqrt(eps*mu);
omega = 2*pi*f;
c = 1/sqrt(mu*eps);
eta = sqrt(mu/eps);

l_m = 1; % Magnetic length of dipoles. This is just a normalization constant.


%% Grab constants and data from CST simulation
CST_data = cst_data_reader('TangentialFieldsFinal.txt');
H_cst = CST_data(1).h_field_center_1_Re + 1i* CST_data(1).h_field_center_1_Im;
H_cst_top = CST_data(3).h_field_top_1_Re + 1i* CST_data(3).h_field_top_1_Im;
waveguide_length = CST_data(1).Par('l_wg') * 1E-3;
waveguide_width = CST_data(1).Par('w_wg') * 1E-3;
waveguide_height = CST_data(1).Par('h_wg') * 1E-3;
waveguide_spacing = CST_data(1).Par('distance_wg') * lambda ;
element_spacing = CST_data(1).Par('distance_slot') * lambda ;
n_elem = 5; % Number of elements per waveguide.

a = waveguide_width;
b = waveguide_height;
l = waveguide_length;

nTx = 2; % number of waveguides
sTx = waveguide_spacing;  % spaceing of waveguides
nDMA = n_elem; % per waveguide (does not have to be the same for all waveguides in general).
sDMA = element_spacing; % DMA element spacing along waveguide

%% Variables and functions
syms x y z X Y Z ... % cartesian coordinates
    m n ... % Mode indexes
    delta_0 dirac ... % a dirac function
    
assume(x,'real'); assume(y,'real'); assume(z,'real');
assume(X,'real'); assume(Y,'real'); assume(Z,'real');
assume(m,{'real','integer'}); assume(n,{'real','integer'});

omega=2*pi*f;

% To insert values, use 'subs(fun, {var1, var2}, {val1, val2})'
%global Rot rot
r = [x;y;z];
ri = [-x;y;-z];
R = [X;Y;Z];
d = norm(r-R,2);
lr =  norm(r,2);
rot = @(V) curl(V,r);
Rot = @(V) curl(V,R);

%% half-space propagation / image theory:
G_air=symZeros(3,3);
G_air(1,1) = 2*exp(-1i*k*sqrt(x^2 + y^2 + z^2))*((((y^2 + z^2)*(x^2 + y^2 + z^2)*k^2)/2 + x^2 - y^2/2 - z^2/2)*sqrt(x^2 + y^2 + z^2) + (x^2 + y^2 + z^2)*(x^2 - y^2/2 - z^2/2)*k*1i)/(2*(x^2 + y^2 + z^2)^3*k^2*pi);

G_airFull=symZeros(3,3);
G_airFull(1,1) =  2*exp(-1i*k*sqrt(x^2 + y^2 + z^2))*(((x^2 + y^2 + z^2)*(y^2 + z^2)*k^2)/2 + sqrt(x^2 + y^2 + z^2)*(2*x^2 - y^2 - z^2)*k*1i/2 - y^2/2 - z^2/2 + x^2)/(2*(x^2 + y^2 + z^2)^(5/2)*pi*k^2);
G_airFull(2,1) =  2*3*(k*sqrt(x^2 + y^2 + z^2)*1i + 1 + (-x^2/3 - y^2/3 - z^2/3)*k^2)*exp(-1i*k*sqrt(x^2 + y^2 + z^2))*x*y/(4*(x^2 + y^2 + z^2)^(5/2)*pi*k^2);
G_airFull(3,1) =  2*(3*exp(-1i*k*sqrt(x^2 + y^2 + z^2))*(k*sqrt(x^2 + y^2 + z^2)*1i + 1 + (-x^2/3 - y^2/3 - z^2/3)*k^2)*x*z)/(4*(x^2 + y^2 + z^2)^(5/2)*pi*k^2);

H_air = matlabFunction(-1i*omega*eps*G_air*[l_m;0;0]); % Function of distance vector in (dx, dy, dz)
H_airFull = matlabFunction(-1i*omega*eps*G_airFull*[l_m;0;0]); % Function of distance vector in (dx, dy, dz)

Y_0 = 2*l_m^2*k*omega*eps/(6*pi);  % Free space above ground real part self-impedance

%% wavenumbers
k_x = m * pi / a; k_y = n * pi / b;
k_c = sqrt(k_x^2 + k_y^2);
k_z = sqrt(k^2 - k_c^2);

% If higher order modes are considered, we must remember to take negative
% imaginary part.
if (orderM > 1) || (orderN > 0)
    k_z = real(sqrt(k^2 - k_c^2)) - 1i*imag(sqrt(k^2 - k_c^2));
end

x_hat = [1; 0; 0]; y_hat = [0; 1; 0]; z_hat = [0; 0; 1];

%% Cartesian vector wave functions

psi_e = @(h) cos(m*pi*x/a)*cos(n*pi*y/b)*exp(1i*h*z);
psi_o = @(h) sin(m*pi*x/a)*sin(n*pi*y/b)*exp(1i*h*z);

N_e = @(h) rot(rot(psi_e(h) * z_hat)) / sqrt(k_c^2 + h^2); N_ex = @(h) subs(N_e(h), {x,y,z},{X,Y,Z});
M_o = @(h) rot(psi_o(h) * z_hat); M_ox = @(h) subs(M_o(h), {x,y,z},{X,Y,Z});

genFunE2Pos = -1i/(a*b) * (2-delta_0)/(k_c^2 * k_z) * ...
    (M_o(-k_z)*M_ox(k_z).' + N_e(-k_z)*N_ex(k_z).');
genFunE2Neg = -1i/(a*b) * (2-delta_0)/(k_c^2 * k_z) * ...
    (M_o(k_z)*M_ox(-k_z).' + N_e(k_z)*N_ex(-k_z).');

pl = exp(-1i*k_z*l)/(2*1i*sin(k_z*l));
plB = exp(1i*k_z*l)/(2*1i*sin(k_z*l));
genFunE2s = -1i/(a*b) * (2-delta_0)/(k_c^2 * k_z) * (...
    M_o(k_z) * (pl*M_ox(k_z).' + pl*M_ox(-k_z).') + ...
    N_e(k_z) * (-pl*N_ex(k_z).' + pl*N_ex(-k_z).') + ...
    M_o(-k_z) * (pl*M_ox(k_z).' + plB*M_ox(-k_z).') + ...
    N_e(-k_z) * (pl*N_ex(k_z).' + (-plB)*N_ex(-k_z).') );

% Here we for simplicity only account for the TE_10 mode.
if orderM == 1 && orderN == 0
    G_e2_positive = simplify(subs(gen_Ge2(orderM,orderN,1, genFunE2Pos, genFunE2Neg, genFunE2s),{x, X},{a/2, a/2}));
    G_e2_negative = simplify(subs(gen_Ge2(orderM,orderN,0, genFunE2Pos, genFunE2Neg, genFunE2s),{x, X},{a/2, a/2}));
else
    G_e2_positive = (subs(gen_Ge2(orderM,orderN,1, genFunE2Pos, genFunE2Neg, genFunE2s),{x, X},{a/2, a/2}));
    G_e2_negative = (subs(gen_Ge2(orderM,orderN,0, genFunE2Pos, genFunE2Neg, genFunE2s),{x, X},{a/2, a/2}));
end
H_wgPos = matlabFunction(-1i*omega*eps*G_e2_positive*[l_m;0;0]);
H_wgNeg = matlabFunction(-1i*omega*eps*G_e2_negative*[l_m;0;0]);



%% Define system
x_dma = (0:(nTx-1)) * sTx + a/2;
y_dma = b;
z_dma = (1:nDMA) * sDMA; z_dma = z_dma - mean(z_dma) + l/2;

[Xt,Yt,Zt] = ndgrid(x_dma, y_dma, z_dma);
xyz_dma = [Xt(:), Yt(:), Zt(:)];

x_tx = (0:(nTx-1)) * sTx + a/2;
y_tx = b/2;
z_tx = 0;
[Xt,Yt,Zt] = ndgrid(x_tx, y_tx, z_tx);
xyz_tx = [Xt(:), Yt(:), Zt(:)];

%% Compute Y_tt
% we assume that there's only a single input per waveguide.
if orderM == 1 && orderN == 0
    Y_tt = eye(size(xyz_tx,1)) * (-[l_m,0,0] * H_wgPos(0,0));
else
    Y_tt = eye(size(xyz_tx,1)) * -[l_m,0,0] * H_wgPos(b/2,0,b/2,0);
end

%% Compute Y_st
Y_st = zeros(size(xyz_dma,1),size(xyz_tx,1));

if orderM == 1 && orderN == 0
    yst = @(z,Z) -[l_m,0,0] * H_wgPos(Z,z);
else
    yst = @(z,Z) -[l_m,0,0] * H_wgPos(b/2,Z,b,z);
end

for idn = 1:size(xyz_dma,1)
    for idm = 1:size(xyz_tx,1)
        if abs(xyz_dma(idn,1) - xyz_tx(idm,1)) <= (a/2)
            Y_st(idn,idm) = yst(xyz_dma(idn,3), xyz_tx(idm,3));
        end
    end
end

%% Compute Y_ss
Y_ss = zeros(size(xyz_dma,1),size(xyz_dma,1));

yssAir = @(x,y,z) -[l_m,0,0] * H_air(x,y,z);
if orderM == 1 && orderN == 0
    yssPos = @(z,Z) -[l_m,0,0] * H_wgPos(Z,z);
    yssNeg = @(z,Z) -[l_m,0,0] * H_wgNeg(Z,z);
else
    yssPos = @(z,Z) -[l_m,0,0] * H_wgPos(b,Z,b,z);
    yssNeg = @(z,Z) -[l_m,0,0] * H_wgNeg(b,Z,b,z);
end

for idn = 1:size(xyz_dma,1) % to
    for idm = 1:size(xyz_dma,1) % from
        if idm == idn
            ent = Y_0 + yssPos(xyz_dma(idn,3),xyz_dma(idn,3));
        elseif abs(xyz_dma(idn,1) - xyz_dma(idm,1)) <= (a/2) % Same waveguide
            if xyz_dma(idn,3) >= xyz_dma(idm,3) % Positive
                ent = yssPos(xyz_dma(idn,3), xyz_dma(idm,3)) + ...
                    yssAir(xyz_dma(idn,1)-xyz_dma(idm,1),xyz_dma(idn,2)-xyz_dma(idm,2),xyz_dma(idn,3)-xyz_dma(idm,3));
            else % Negative
                ent = yssNeg(xyz_dma(idn,3), xyz_dma(idm,3)) + ...
                    yssAir(xyz_dma(idn,1)-xyz_dma(idm,1),xyz_dma(idn,2)-xyz_dma(idm,2),xyz_dma(idn,3)-xyz_dma(idm,3));
            end
        else % other waveguide
            ent = yssAir(xyz_dma(idn,1)-xyz_dma(idm,1),xyz_dma(idn,2)-xyz_dma(idm,2),xyz_dma(idn,3)-xyz_dma(idm,3));
        end
        Y_ss(idn,idm) = ent;
    end
end

%% Compute currents and fields:
zm = (CST_data(1).Z_mm -CST_data(1).Z_mm(1)) * 1E-3; % Measurement points

Y_s0 = eye(size(xyz_dma,1))*(2-15.7934i) * l_m^2;

Yp = Y_tt- Y_st.' * inv(Y_ss+Y_s0) * Y_st;
Y0 = ones(nTx,1) * Y_intrinsic_source*l_m^2;

j = [1;1];
Yin = sum(Yp * j,2) ./ j;
Gamma = diag( -(Yin - Y0) ./ (Yin + Y0) );
j_t = (eye(nTx) + Gamma) * j;

% normalize power:
Ps = real(j_t' * inv((eye(nTx) - Gamma' * Gamma)) * Yp * j_t) / 2;

j = j / sqrt(Ps);
j_t = (eye(nTx)+Gamma) * j ;
j_s = -(Y_ss+Y_s0)^(-1)* Y_st * j_t; % DMA element currents

Ps = real(j_t' * inv((eye(nTx) - Gamma' * Gamma)) * Yp * j_t) / 2; % Power supplied
Pt = real(j_t' * Yp * j_t) / 2;  % Power radiated

if plot_stuff
    H = zeros(length(zm),1);
    H_top = zeros(length(zm),1);
    x0 = a/2;  % The center of the waveguide of interest
    for n = 1:length(zm)
        
        H_field = zeros(3,1);
        H_field_top = zeros(3,1);
        
        for id_tx = 1:size(xyz_tx,1)
            if xyz_tx(id_tx,1) == x0 % If it's in the waveguide of interest
                H_field = H_field + j_t(id_tx) * ...
                    (H_wgPos(xyz_tx(id_tx,3),zm(n))*(zm(n)>=xyz_tx(id_tx,3)) + H_wgNeg(xyz_tx(id_tx,3),zm(n))*((zm(n)<xyz_tx(id_tx,3))));
            end
        end
        for id_dma = 1:size(xyz_dma,1)
            if xyz_dma(id_dma,1) == x0 % If it's in the waveguide of interest
                H_field = H_field + j_s(id_dma) * ...
                    (H_wgPos(xyz_dma(id_dma,3),zm(n))*(zm(n)>=xyz_dma(id_dma,3)) + H_wgNeg(xyz_dma(id_dma,3),zm(n))*((zm(n)<xyz_dma(id_dma,3))));
            end
            
            H_field_top = H_field_top - j_s(id_dma) * H_air(x0 - xyz_dma(id_dma,1), 0, zm(n)-xyz_dma(id_dma,3));
        end
        
        H(n) = H_field(1);
        H_top(n) = H_field_top(1);
        
    end

    fig = figure; hold on;
    
    L(1) = plot(nan, nan, 'k-');%'LineWidth',1
    L(2) = plot(nan, nan, 'k:*');
    legend(L, {'CST', 'Model'})
    
    zaxis = zm * 1E3;
    yyaxis left; hold on; grid on; l1=legend;
    ylabel('$|(\mathbf{h})_3|$  $[\mathrm{A}\,\mathrm{m}^{-1}]$')
    xlabel('$x$ $[\mathrm{mm}]$')
    plot(zaxis,abs(H_cst),'color',colours(1,:),'displayName','CST','MarkerIndices',unique(round(linspace(1,length(zaxis),25),0)))
    plot(zaxis,abs(H),':*','color',colours(1,:),'displayName','Model','MarkerIndices',unique(round(linspace(1,length(zaxis),20),0)))
    xlim([min(zaxis),max(zaxis)])
    
    yyaxis right; hold on; grid on;
    ylabel('$\angle (\mathbf{h})_3$  $[\mathrm{rad}]$')
    xlabel('$x$ $[\mathrm{mm}]$')
    plot(zaxis,(angle(H_cst)),'color',colours(2,:),'displayName','CST','MarkerIndices',unique(round(linspace(1,length(zaxis),25),0)))
    plot(zaxis,(angle(H)),':*','color',colours(2,:),'displayName','Model','MarkerIndices',unique(round(linspace(1,length(zaxis),20),0)))
    xlim([min(zaxis),max(zaxis)])
    
    ax = gca;
    ax.YAxis(1).Color = colours(1,:);
    ax.YAxis(2).Color = colours(2,:);
    leg = ax.Legend;
    leg.String = leg.String(1:2);
    
    set(0, 'CurrentFigure', fig)
    
end



%% Farfield pattern

farfield_data = cst_farfield_reader('farfieldPatternRotated.txt');

distance = 1E10;
phi = deg2rad(linspace(0,180,180)); % Only the half-space for positive y
theta = linspace(0,pi,100);
[P,T] = meshgrid(phi,theta);
G = zeros(size(phi,2), size(theta,2));

for id_ph = 1:length(phi)
    ph = phi(id_ph);
    for id_th = 1:length(theta)
        th = theta(id_th);
        H = zeros(3,1);
        for id_dma = 1:1:size(xyz_dma,1)
            r = [cos(th); sin(ph)*sin(th); cos(ph)*sin(th)] * distance;
            H = H - j_s(id_dma) * H_airFull(r(1) - xyz_dma(id_dma,1), r(2)-0, r(3)-xyz_dma(id_dma,3));
        end
        
        TransformSpherical = [sin(th)*cos(ph), sin(th)*sin(ph), cos(th);...
            cos(th)*cos(ph), cos(th)*sin(ph), -sin(th); ...
            -sin(ph), cos(ph), 0];
        TransformSpherical = [TransformSpherical(:,1), TransformSpherical(:,2), TransformSpherical(:,3)];
        
        H_sph = TransformSpherical * H;
        
        W = eta * (H_sph' * H_sph) /2;
        U = distance^2 * W;
        G(id_ph,id_th) = (4*pi/Ps) * U; % U = radiation intensity W/unit solid angle
    end
end

GdB = 10*log10(G);

figureScale = 2; % Scales the figure size and font size.
figureHeight = 12;
figure('units','centimeters','position',[0,0,linewidth*figureScale, figureHeight]);

titleHeight = 0.93; % Moves title above radiation patterns
scalingMinimumRadius = 0;
xyScale = 1; % Putting this higher than 1 might cause figure to clip out of plot

scaleMax = round(max([GdB(:);farfield_data.grlz_db]));
scaleMin = -30;
dynRange = scaleMax - scaleMin;


%% Plot CST simulated pattern:
G1 = farfield_data.grlz_db;
theta = farfield_data.theta;
phi = farfield_data.phi;

uTh = unique(theta);
uPh = unique(phi);

hax = subplot(1,2,1); hold on;
hplot = patternCustom(G1,rad2deg(theta),rad2deg(phi),'CoordinateSystem','polar'); % CoordinateSystem = polar | rectangular,
view(160,30);
titleStr = ['$G_\mathrm{CST}$ [dBi]'];
t = text(0.5,titleHeight,titleStr,'interpreter','latex','units','normalized','HorizontalAlignment','center');
ca = gca;
ca.CameraViewAngleMode = 'manual';
ca.DataAspectRatioMode = 'manual';
ca.PlotBoxAspectRatioMode = 'manual';
ca1 = ca;
c = colorbar;
for n = 1:length(uTh)
    th = uTh(n);
    Gs = (G1(theta == th) - scaleMin) / dynRange;
    Gs(Gs < scalingMinimumRadius) = scalingMinimumRadius;
    
    hplot.XData(:,n) = cos(uPh).*sin(th) .* Gs * xyScale;
    hplot.YData(:,n) = sin(uPh).*sin(th) .* Gs * xyScale;
    hplot.ZData(:,n) = cos(th) .* Gs;
end
caxis([scaleMin, scaleMax])
ldiff = c.Position(4) - 0.6;
c.Position(4) = c.Position(4) - ldiff;
c.Position(2) = c.Position(2) + ldiff/2;
c.Position(1) = c.Position(1) + 0.03;
c.Visible = 'off';
c1 = c;
removeGuides(hplot, remove_spherical_guides, remove_cartesian_guides)
hax.CameraPosition = [6.7831,   15.7425,    6.9781];
hax.DataAspectRatio = [1 1 1];
hax.PlotBoxAspectRatio = [1,1,1];
xlim([-1.2000 1.2000]); zlim([-1.2000 1.2000]); ylim([-1.2000 1.2000]);



%% Plot Model simulated pattern:
G2 = GdB'; G2 = G2(:);
theta = T(:);
phi = P(:);

uTh = unique(theta);
uPh = unique(phi);

hax = subplot(1,2,2); hold on;
hplot = patternCustom(G2,rad2deg(theta),rad2deg(phi),'CoordinateSystem','polar'); % CoordinateSystem = polar | rectangular,
view(160,30);
titleStr = ['$G_\mathrm{model}$ [dBi]'];
t = text(0.5,titleHeight,titleStr,'interpreter','latex','units','normalized','HorizontalAlignment','center');
ca = gca;
ca.CameraViewAngleMode = 'manual';
ca.DataAspectRatioMode = 'manual';
ca.PlotBoxAspectRatioMode = 'manual';
ca2 = ca;
c = colorbar;
for n = 1:length(uTh)
    th = uTh(n);
    Gs = (G2(theta == th) - scaleMin) / dynRange;
    Gs(Gs < scalingMinimumRadius) = scalingMinimumRadius;
    
    hplot.XData(:,n) = cos(uPh).*sin(th) .* Gs * xyScale;
    hplot.YData(:,n) = sin(uPh).*sin(th) .* Gs * xyScale;
    hplot.ZData(:,n) = cos(th) .* Gs;
end
caxis([scaleMin, scaleMax]);

ldiff = c.Position(4) - 0.6;
c.Position(4) = c.Position(4) - ldiff;
c.Position(2) = c.Position(2) + ldiff/2;
c.Position(1) = c.Position(1) ;
c2 = c;
removeGuides(hplot, remove_spherical_guides, remove_cartesian_guides)
hax.CameraPosition = [6.7831,   15.7425,    6.9781];
hax.DataAspectRatio = [1 1 1];
hax.PlotBoxAspectRatio = [1,1,1];
xlim([-1.2000 1.2000]); zlim([-1.2000 1.2000]); ylim([-1.2000 1.2000]);
drawnow;

curPos = ca2.Position;
curPos(1) = curPos(1) - 0.06;
ca2.Position = curPos;

curPos = ca1.Position;
curPos(1) = curPos(1) - 0.01;
ca1.Position = curPos;


%% Plot farfield cuts.

figure; hold on; grid; l = legend('location','southoutside');
l.NumColumns = 2;

phi = P(1,:)';
theta = T(:,1);
scanAngle = pi/2;
[~,id_th] = min(abs(theta - scanAngle));

colourIndx = colourIndx+1; styleIndx = colourIndx;
lineStyle = [lineSpec{mod(styleIndx,length(lineSpec))+1},markerSpec{mod(styleIndx,length(markerSpec))+1}];
plot(phi,GdB(:,id_th),lineStyle,'color',colours(colourIndx,:),'displayName','Model, $\theta = \frac{\pi}{2}\,\mathrm{rad}$','MarkerIndices',unique(round(linspace(1,length(phi),20),0)))

tf = round(farfield_data.theta,5) == round(scanAngle,5);
colourIndx = colourIndx+1; styleIndx = colourIndx;
lineStyle = [lineSpec{mod(styleIndx,length(lineSpec))+1},markerSpec{mod(styleIndx,length(markerSpec))+1}];
plot(farfield_data.phi(tf),farfield_data.grlz_db(tf),lineStyle,'color',colours(colourIndx,:),'displayName','CST, $\theta = \frac{\pi}{2}\,\mathrm{rad}$','MarkerIndices',unique(round(linspace(1,length(farfield_data.phi(tf)),20),0)))

scanAngle = deg2rad(45);
[m,id_th] = min(abs(theta - scanAngle));
colourIndx = colourIndx+1; styleIndx = colourIndx;
lineStyle = [lineSpec{mod(styleIndx,length(lineSpec))+1},markerSpec{mod(styleIndx,length(markerSpec))+1}];
plot(phi,GdB(:,id_th),lineStyle,'color',colours(colourIndx,:),'displayName','Model, $\theta = \frac{\pi}{4}\,\mathrm{rad}$','MarkerIndices',unique(round(linspace(1,length(phi),20),0)))

tf = round(farfield_data.theta,5) == round(scanAngle,5);
colourIndx = colourIndx+1; styleIndx = colourIndx;
lineStyle = [lineSpec{mod(styleIndx,length(lineSpec))+1},markerSpec{mod(styleIndx,length(markerSpec))+1}];
plot(farfield_data.phi(tf),farfield_data.grlz_db(tf),lineStyle,'color',colours(colourIndx,:),'displayName','CST, $\theta = \frac{\pi}{4}\,\mathrm{rad}$','MarkerIndices',unique(round(linspace(1,length(farfield_data.phi(tf)),20),0)))

ylabel('Gain [dBi]');
xlabel('$\phi$ [rad]')
ylim([-20,15])
xlim([phi(2), 3.12414])

%% Plot amplitude and phase-dependencies

c = linspace(-5,5,100);
vartheta = 1./(1 + 1i*c);

fig = figure; hold on;

yyaxis left; hold on; grid on;
ylabel('$\|\vartheta\|$  $[\Omega]$')
xlabel('$c$ $[\mathrm{S}]$')
plot(c,abs(vartheta),'color',colours(1,:),'displayName','CST','MarkerIndices',unique(round(linspace(1,length(c),25),0)))
xlim([min(c),max(c)])
xticks(-5:5)
ylim([0, 1]);
yticks(0:0.25:1)

yyaxis right; hold on; grid on;
ylabel('$\angle\vartheta$  $[\mathrm{rad}]$')
xlabel('$c$ $[\mathrm{S}]$')
plot(c,angle(vartheta),'color',colours(2,:),'displayName','CST','MarkerIndices',unique(round(linspace(1,length(c),25),0)))
xlim([min(c),max(c)])
ylim([-0.5*pi, 0.5*pi]);
yticks([-0.5*pi -0.25*pi 0 0.25*pi 0.5*pi])
yticklabels({'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'})

ax = gca;
ax.YAxis(1).Color = colours(1,:);
ax.YAxis(2).Color = colours(2,:);

set(0, 'CurrentFigure', fig)

%% Functions
function G = gen_Ge2(mm,nn,positive, genFunE2Pos, genFunE2Neg, genFunE2s)
syms z Z k dirac m n delta_0
z_hat = [0; 0; 1];

G =-1/(k^2)*dirac * (z_hat * z_hat');
for o = 0:nn
    for p = 0:mm
        if ~(o==0 && p == 0)
            d_0 = (o == 0) || (p == 0);
            if positive
                G = G + subs(genFunE2Pos + genFunE2s, {m,n,delta_0},{p,o,d_0});
            else
                G = G + subs(genFunE2Neg + genFunE2s, {m,n,delta_0},{p,o,d_0});
            end
        end
    end
end
end

function A = symZeros(n,m)
A = sym('A', [n,m]);
for idn = 1:n
    for idm = 1:m
        A(idn,idm) = 0;
    end
end
end

function removeGuides(hplot, remove_spherical_guides, remove_cartesian_guides)
if remove_spherical_guides || remove_cartesian_guides
    for n = 1:length(hplot.Parent.Children)
        if remove_cartesian_guides
            if strcmp(hplot.Parent.Children(n).Type,'line')
                if length(hplot.Parent.Children(n).XData) == 2
                    hplot.Parent.Children(n).Visible = 'off';
                end
            elseif strcmp(hplot.Parent.Children(n).Type,'text')
                if length(hplot.Parent.Children(n).String) == 1
                    hplot.Parent.Children(n).Visible = 'off';
                end
            end
        end
        
        if remove_spherical_guides
            if strcmp(hplot.Parent.Children(n).Type,'line')
                if length(hplot.Parent.Children(n).XData) ~= 2
                    hplot.Parent.Children(n).Visible = 'off';
                end
            elseif strcmp(hplot.Parent.Children(n).Type,'text')
                if length(hplot.Parent.Children(n).String) == 4
                    hplot.Parent.Children(n).Visible = 'off';
                end
            elseif strcmp(hplot.Parent.Children(n).Type,'patch')
                hplot.Parent.Children(n).Visible = 'off';
            end
        end
    end
end
end

function setFigOptions()
%CHANGEFIDEFAULTS Call before creating figures
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontSize',10)
set(0,'DefaultAxesFontName', 'Latin Modern Roman')
set(0,'DefaultTextFontname', 'Latin Modern Roman')
end

function data = cst_data_reader(file)
%Imporing CST exports to matlab.

path = '';
if exist(file, 'file')
    if ~strcmpi(file(end-3:end),'.txt')
        file = [file,'.txt'];
    end
end
if ~exist(file, 'file')
    [file,path] = uigetfile({'*.txt','CST PLOTDATA Export (*.txt)'});
end

data = [];
dataIndx = 0;
saveIndx = {};

if exist(strrep(file,'txt','mat'), 'file')
    load(strrep(file,'txt','mat'))
else
    fid = fopen([path,file]);
    tab = sprintf('\t');
    paramExp = '(?<var>[a-z_A-Z0-9]+?)=(?<val>-?[0-9.]*)[; ]*';
    nameExp = '(?<var>[a-zA-Z0-9_,. \-]+) [/\[\(]';
    formatExp = '((\[)|(/ ))(?<format>[a-zA-Z.]+)';
    
    tline = fgetl(fid); % Read first line
    while (tline ~= -1) % Check if there is more data
        if tline(1) == '#' % This is new header
            dataIndx = dataIndx + 1;
            
            % Realdall parameters:
            out = regexp(tline,paramExp,'tokens');
            data(dataIndx).Par = containers.Map;
            if ~isempty(out)
                for parNum = 1:length(out)
                    data(dataIndx).Par(out{parNum}{1}) = str2double(out{parNum}{2});
                end
            end
            
            % Create arrays for all columns.
            tline = fgetl(fid);
            saveIndx = {}; % Contains data column to array name link
            if tline == -1
                error('File ended unexpectingly, no data in file.');
            end
            C = strsplit(tline(2:end),tab);
            for CIndx = 1:length(C)
                outName = regexp(C{CIndx},nameExp,'tokens');
                outFormat = regexp(C{CIndx},formatExp,'tokens');
                name = outName{1}{1};
                if name(end) == ' '
                    name(end) = [];
                end
                if name(end) == '.'
                    name(end) = [];
                end
                if isempty(outFormat)
                    format = [];
                else
                    format = outFormat{1}{2};
                end
                
                name = strrep(name,',','_');
                name = strrep(name,' ','_');
                name = strrep(name,'.','_');
                name = strrep(name,'-','_');
                format = strrep(format,',','_');
                format = strrep(format,' ','_');
                format = strrep(format,'.','_');
                format = strrep(format,'-','_');
                name = [name , '_', format];
                saveIndx{1,end+1} = name;
                saveIndx{2,end} = CIndx;
                eval(['data(dataIndx).', name,' = [];']); % Create variable to store data
            end
            
            fgetl(fid); % Third line of header is discarded.
        else % Read data
            C = strsplit(tline(1:end),tab);
            for SI_indx = 1:size(saveIndx,2)
                eval(['data(dataIndx).', saveIndx{1,SI_indx}, ...
                    ' = [data(dataIndx).', saveIndx{1,SI_indx},'; ', C{saveIndx{2,SI_indx}},'];']);
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    save(strrep(file,'txt','mat'),'data')
end
end

function [data] = cst_farfield_reader(file)
%Imporing CST farfield plots

if exist(file, 'file')
    if ~strcmpi(file(end-3:end),'.txt')
        file = [file,'.txt'];
    end
    rad_pat = importdata(file);
else
    [file,path] = uigetfile({'*.txt','CST FARFIELD Export (*.txt)'});
    rad_pat = importdata([path,file]);
end

theta = deg2rad(rad_pat.data(:,1));
phi = deg2rad(rad_pat.data(:,2));
grlz_db = rad_pat.data(:,3);
%grlz_theta_db = rad_pat.data(:,4);
%phase_theta = deg2rad(rad_pat.data(:,5));
%grlz_phi_db = rad_pat.data(:,6);
%phase_phi = deg2rad(rad_pat.data(:,7));
axRatio = rad_pat.data(:,8);

data.theta =   [theta; theta(phi == 0)];
data.phi =     [phi; deg2rad(360) + phi(phi == 0)];
data.grlz_db = [grlz_db; grlz_db(phi == 0)];
data.grlz_lin = 10.^(data.grlz_db/10);
data.axRatio = [axRatio; axRatio(phi == 0)];
end
