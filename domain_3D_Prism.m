%% Function is fed with inputs from MoL_Script
% Note: This function is fed with dimensionless paramters!!!

function [fval] = domain_3D_Prism(t,y0,X0,Y0,Z0,P1,P2_x,P2_y,P2_z,P3)
   
    % Define the number of spatial 'slices'  
    X = X0+1; % step size = 1/X in X-direction
    Y = Y0+1; % step size = 1/Y in Y-direction
    Z = Z0+1; % step size = 1/Z in Z-direction
    % Restructure input (y) back to matrix
    sol = reshape(y0,X0+Y0+Z0,Y0,Z0); 

    % Define the [Drug] vector
    Drug = zeros(X+1,Y+1,Z+1); % Preallocate space in a vector
    % Extract [Drug] values from input
        % Rows = Z-direction (N0 slices)
        % Cols = R-direction (R0 slices)
    Drug(2:X,2:Y,2:Z) = sol(1:X0,1:end,1:end); 
    % d[Drug]_dx = 0 at x = 0 
    Drug(1,1:end,1:end) = (1/3)*(-Drug(3,1:end,1:end)+...
        4*Drug(2,1:end,1:end)); 
    % d[Drug]_dy = 0 at y = 0 
    Drug(1:end,1,1:end) = (1/3)*(-Drug(1:end,3,1:end)+...
        4*Drug(1:end,2,1:end));
    % d[Drug]_dz = 0 at z = 0 
    Drug(1:end,1:end,1) = (1/3)*(-Drug(1:end,1:end,3)+...
        4*Drug(1:end,1:end,2));
    % [Drug] outside hydrogel is zero
    Drug(X+1,1:end,1:end) = 0; 
    Drug(1:end,Y+1,1:end) = 0;
    Drug(1:end,1:end,Z+1) = 0;
    
    % Define the [Drug-CD] vector
    Drug_CD = zeros(X+1,Y+1,Z+1); % Preallocate space in a vector
    % Extract [Drug-CD] values from input
        % Rows = Z-direction (N0 slices)
        % Cols = R-direction (R0 slices)
    Drug_CD(2:X,2:Y,2:Z) = sol(X0+1:2*X0,1:end,1:end);
    % d[Drug-CD]_dx = 0 at x = 0 
    Drug_CD(1,1:end,1:end) = (1/3)*(-Drug_CD(3,1:end,1:end)+...
        4*Drug_CD(2,1:end,1:end)); 
    % d[Drug-CD]_dy = 0 at y = 0 
    Drug_CD(1:end,1,1:end) = (1/3)*(-Drug_CD(1:end,3,1:end)+...
        4*Drug_CD(1:end,2,1:end));
    % d[Drug-CD]_dz = 0 at z = 0 
    Drug_CD(1:end,1:end,1) = (1/3)*(-Drug_CD(1:end,1:end,3)+...
        4*Drug_CD(1:end,1:end,2));
    % [Drug-CD] outside hydrogel is zero
    Drug_CD(X+1,:,:) = 0;
    Drug_CD(:,Y+1,:) = 0;
    Drug_CD(:,:,Z+1) = 0;

    % Define the derivatives
    d2Drug_dx2 = zeros(size(Drug));  % Preallocate space in vectors
    d2Drug_dy2 = zeros(size(Drug)); 
    d2Drug_dz2 = zeros(size(Drug)); 
    dDrugCD_dt = zeros(size(Drug_CD));
    dC_dt = zeros(size(Drug));
    dMR_dt = zeros(X0,Y0,Z0);
   % dMR_dt = -1/2*P2*(Drug(N+1,)-Drug(N))/(1/N);

    for i = 2:X  % Iterate through X-direciton
        for j = 2:Y  % Iterate through Y-direction
            for k = 2:Z % Iterate through Z-direction
                
                % Rb,i* = Net Rate of Drug-CD Binding 
                dDrugCD_dt(i,j,k) = P1*Drug(i,j,k)*(P3-Drug_CD(i,j,k))-...
                    Drug_CD(i,j,k);
                
                if i > 2 && i < X  % If slice is inbetween edge and center 
                    d2Drug_dx2(i,j,k) = (Drug(i+1,j,k)-2*Drug(i,j,k)+...
                        Drug(i-1,j,k))/((1/X0)^2); % Change in X-direction
                elseif i == 2 % If slice is on center 
                     d2Drug_dx2(i,j,k) = (Drug(i,j,k)-Drug(i-1,j,k))/...
                        ((1/X0)^2);     % Change in X-direction
                else % i == X % If slice is on edge 
                    d2Drug_dx2(i,j,k) = (-2*Drug(i,j,k)+Drug(i-1,j,k))/...
                        ((1/X0)^2);     % Change in X-direction
                end 
                
                if j > 2 && j < Y % If slice is inbetween edge and center 
                    d2Drug_dy2(i,j,k) = (Drug(i,j+1,k)-2*Drug(i,j,k)+...
                        Drug(i,j-1,k))/((1/Y0)^2); % Change in Y-direction
                elseif j == 2 % If slice is on center 
                    d2Drug_dy2(i,j,k) = (Drug(i,j,k)-Drug(i,j-1,k))/...
                        ((1/Y0)^2); % Change in Y-direction
                else  % j == Y % If slice is on edge 
                    d2Drug_dy2(i,j,k) = (-2*Drug(i,j,k)+Drug(i,j-1,k))/...
                        ((1/Y0)^2); % Change in Y-direction
                end 
                    
                if k > 2 && k < Z % If slice is inbetween edge and center 
                    d2Drug_dz2(i,j,k) = (Drug(i,j,k+1)-2*Drug(i,j,k)+...
                        Drug(i,j,k-1))/((1/Z0)^2); % Change in Z-direction
                elseif k == 2 % If slice is on center 
                    d2Drug_dz2(i,j,k) = (Drug(i,j,k)-Drug(i,j,k-1))/...
                        ((1/Z0)^2); % Change in Z-direction
                else % k == Z % If slice is on edge 
                    d2Drug_dz2(i,j,k) = (-2*Drug(i,j,k)+Drug(i,j,k-1))/...
                        ((1/Z0)^2); % Change in Z-direction
                end 
                
                % Change in [Drug] in the hydrogel at spatial slice (i) 
                    dC_dt(i,j,k) = P2_x*d2Drug_dx2(i,j,k) + ...
                        P2_y*d2Drug_dy2(i,j,k) + P2_z*d2Drug_dz2(i,j,k)...
                        - dDrugCD_dt(i,j,k);
               
            end 
        end
    end 
    
    % Extract fval as a vector of all other vectors
    % Note: The format of fval should be identical to the input (y)
    % Calculate change in 'Moles Released' in x-dir
    dMR_dt(1,:,:) = -1/2*P2_x*(Drug(X+1,2:Y,2:Z)-Drug(X,2:Y,2:Z))/(1/X0);
    % Calculate change in 'Moles Released' in y-dir
    dMR_dt(:,1,:) = -1/2*P2_y*(Drug(2:X,Y+1,2:Z)-Drug(2:X,Y,2:Z))/(1/Y0);
    % Calculate change in 'Moles Released' in z-dir
    dMR_dt(:,:,1) = -1/2*P2_z*(Drug(2:X,2:Y,Z+1)-Drug(2:X,2:Y,Z))/(1/Z0);
    % Output matrix: [ d[Drug]/dt; d[Drug-CD]/dt; d('Moles Released')/dt
    fval = [dC_dt(2:X,2:Y,2:Z); dDrugCD_dt(2:X,2:Y,2:Z); dMR_dt];
    fval = reshape(fval,numel(fval),1); % Restructure output as vector
end
