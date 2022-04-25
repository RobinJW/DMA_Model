function [Y_rr] = Coupling_Dipoles(f, l, xyz_user, mu, epsilon)

%% Physical constants - Dont edit
k = 2*pi*f*sqrt(epsilon*mu);

%% Y_rr calculation 
Ga_e2 = @(r,rhat) ((norm(r-rhat)^2 - (r(3)-rhat(3))^2)/norm(r-rhat)^2 - ...
                    1i*(norm(r-rhat)^2 - 3*(r(3)-rhat(3))^2)/(norm(r-rhat)^3*k)...
                    - (norm(r-rhat)^2 - 3*(r(3)-rhat(3))^2)/(norm(r-rhat)^4*k^2))...
                    *exp(-1i*k*norm(r-rhat))/(4*pi*norm(r-rhat));
                
Y_rr = zeros(size(xyz_user,1),size(xyz_user,1));

for row = 1:size(xyz_user,1)
   for col = 1:size(xyz_user,1)
      % If diagonal element -> self-admittance
      if(row == col)
          Y_rr(row,col) = l^2*k*2*pi*f*epsilon/(6*pi);
      else 
          Y_rr(row,col) = l^2*1i*2*pi*f*epsilon*...
              Ga_e2(xyz_user(row,:),xyz_user(col,:));
      end
   end
end

end

