function [Y_rs] = GenChannel(channel_type, lambda, ant_xyz, user_xyz)

    L = size(ant_xyz,1);                % Number of DMA elements
    M = size(user_xyz,1);               % Number of users
    
    k = 2*pi/lambda;

    % Rayleigh fading 
    if channel_type  
        
        % Distances between users and DMA elements
        D = zeros(M,L);
        for idn = 1:L
            for idm = 1:M
                D(idm,idn) = norm(ant_xyz(idn,:) - user_xyz(idm,:),2);
            end
        end
        % Pathloss
        PL = (lambda ./ (4*pi*D)).^2;
        
        % Uncorrelated Complex Gaussian realizations
        Yrs_uncorr = (randn(M, L) + 1i * randn(M, L)) .* sqrt(PL/2);

        % Compute correlation coefficient.
        dZ = squareform(pdist(ant_xyz(:,3)));
        dR = squareform(pdist(ant_xyz));
        Sigma = 3./2.*((1 + (-k.^2.*dZ.^2 - 1)./(dR.^2.*k.^2) + ...
                 3.*dZ.^2./(dR.^4.*k.^2)).*sin(k.*dR)./(k.*dR) +...
                 cos(k.*dR).*(1./(k.*dR) - 3.*dZ.^2./(dR.^3.*k))./(k.*dR)); 

        Sigma(isnan(Sigma)) = 1;
        sq_Sigma = real(Sigma^(1/2)); % Real operator due to imaginary part
                                      % being a product of quantization

        Y_rs = Yrs_uncorr * sq_Sigma;
        
    % LoS channel    
    else            

        % Distances between users and DMA elements
        D = zeros(M, L);
        dz = zeros(M, L);
        for idn = 1:L
            for idm = 1:M
                D(idm,idn) = norm(ant_xyz(idn,:) - user_xyz(idm,:),2);
                dz(idm,idn) = abs(ant_xyz(idn,3) - ant_xyz(idm,3));
            end
        end
        % Polar angle
        theta = pi/2 - asin(dz./D);

        % Pathloss
        PL = (lambda ./ (4*pi*D)).^2.* (3/2 * sin(theta).^2) .* ...
            (6/2 * sin(theta).^2);
        Y_rs = sqrt(PL) .* exp(-1i*k*D);

    end

end

