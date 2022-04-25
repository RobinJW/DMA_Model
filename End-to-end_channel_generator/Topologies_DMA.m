function [ant_xyz, rf_xyz] = Topologies_DMA(site_xyz,...
                     N, Lmu, wvg_spacing, elem_spacing, S_mu, a, b, Plot_topology)
        
    % Total number of antennas
    L = Lmu*N;
    
    % Pre-allocating
    ant_xyz = zeros(L*size(site_xyz,1),3);
    rf_xyz = zeros(N*size(site_xyz,1),3);
    
    % The sites are generated according to a local coordinate system, and
    % then shifted according to site_xyz
    for ksite = 1:size(site_xyz,1)
        % Coordinates of RF chains
        z_rf = (0:N-1)*wvg_spacing + a/2; 
        y_rf = b/2;
        x_rf = 0;

        % Coordinates of DMA elements
        z_dma = z_rf;
        y_dma = b;
        x_dma = (1:Lmu) * elem_spacing;
        x_dma = x_dma - mean(x_dma) + S_mu/2;
        
        % Store coordinates for antennas in site ksite
        [Xant,Yant,Zant] = ndgrid(x_dma+site_xyz(ksite,1),...
                               y_dma+site_xyz(ksite,2),...
                               z_dma+site_xyz(ksite,3));
                           
        ant_xyz((L*(ksite-1))+1:(L*ksite),:) = ...
                [Xant(:), Yant(:), Zant(:)];
            
        % Store coordinates for RF chains in site ksite
        [Xrf,Yrf,Zrf] = ndgrid(x_rf+site_xyz(ksite,1),...
                       y_rf+site_xyz(ksite,2),...
                       z_rf+site_xyz(ksite,3));
        rf_xyz((N*(ksite-1))+1:(N*ksite),:) = ...
                [Xrf(:), Yrf(:), Zrf(:)];
                
    end
    
    %-----------------------------------------------------------------
    % Plotting deployment
    %-----------------------------------------------------------------
    if (~Plot_topology)
        return
    end
    
    figure(); hold on; grid on; 
    xlabel('x'); ylabel('y'); zlabel('z');
    antplt = scatter3(ant_xyz(:,1),ant_xyz(:,2),ant_xyz(:,3));
    rfplt = scatter3(rf_xyz(:,1),rf_xyz(:,2),rf_xyz(:,3));
    
    % Loop to plot the waveguides
    xtot = zeros(4*N,4);
    ytot = zeros(0,4);
    ztot = zeros(0,4);
    xyz_aux = rf_xyz;

    for idn = 1:N
        zt = [0, a, a, 0; ...
            0, a, a, 0; ...
            0, 0, 0, 0; ...
            a, a, a ,a] + (xyz_aux(2,3)-xyz_aux(1,3))*(idn-1);
        yt = [0, 0, 0, 0; ...
            b, b, b, b; ...
            0, b, b, 0; ...
            0, b, b, 0];
        xt = [0, 0, S_mu, S_mu; ...
            0, 0, S_mu, S_mu;...
            0, 0, S_mu, S_mu;...
            0, 0, S_mu, S_mu] + (xyz_aux(2,1)-xyz_aux(1,1))*(idn-1);

        xtot(4*(idn-1)+1 : 4*idn, :) = xt + site_xyz(1,1);
        ytot(4*(idn-1)+1 : 4*idn, :) = yt + site_xyz(1,2);
        ztot(4*(idn-1)+1 : 4*idn, :) = zt + site_xyz(1,3);

    end
    h1 = patch('XData',xtot','YData',ytot','ZData',ztot');
    h1.FaceAlpha = 0.2;

    axis([0 2*S_mu 0 0.2 site_xyz(1,3) site_xyz(1,3)+(N+1)*wvg_spacing])
    legend([antplt, rfplt, h1],{'DMA elements',  'RF chains', 'Waveguides'})
    title('Deployment DMA')
    view([135,37]);
end