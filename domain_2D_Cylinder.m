%% Function is fed with inputs from MoL_Script
% Note: This function is fed with dimensionless paramters!!!

function [fval] = domain_2D_Cylinder(t,y0,N0,R0,P1,P2_z,P2_r,P3)
   
    % Define the number of spatial 'slices'  
    N = N0+1; % step size = 1/N in z-direction
    R = R0+1; % step size = 1/R in r-direction
    sol = reshape(y0,N0*2+2,N0);    % Restructure input (y) back to matrix

    % Define the [Drug] vector
    Drug = zeros(N+1,R+1); % Preallocate space in a vector
    Drug_CD = Drug;        % Preallocate space in a vector
    Drug(2:N,2:R) = sol(1:N0,1:end); % Extract [Drug] values from input
                          % Rows = Z-direction (N0 slices)
                          % Cols = R-direction (R0 slices)
    Drug(1,1:end) = (1/3)*(-Drug(3,:)+4*Drug(2,:)); % dCL_dz = 0 at z = 0
    Drug(1:end,1) = (1/3)*(-Drug(:,3)+4*Drug(:,2)); % dCL_dr = 0 at r = 0
    % [Drug] outside hydrogel is zero
    Drug(N+1,1:end) = 0;
    Drug(:,R+1) = 0;
  
    % Define the [Drug-CD] vector
    Drug_CD(2:N,2:R) = sol(N0+1:2*N0,:);% Extract [Drug] values from input
                          % Rows = Z-direction (N0 slices)
                          % Cols = R-direction (R0 slices)
    % dCL_dz = 0 at z = 0
    Drug_CD(1,1:end) = (1/3)*(-Drug_CD(3,:)+4*Drug_CD(2,:));
    % dCL_dr = 0 at r = 0
    Drug_CD(1:end,1) = (1/3)*(-Drug_CD(:,3)+4*Drug_CD(:,2));
    % [Drug-CD] outside hydrogel is zero
    Drug_CD(N+1,1:end) = 0;     
    Drug_CD(1:end,R+1) = 0;

    % Define the derivatives
    d2Drug_dz2 = zeros(size(Drug));  % Preallocate space in vectors
    d2Drug_dr2 = zeros(size(Drug)); 
    dDrugCD_dt = zeros(size(Drug_CD));
    dC_dt = zeros(size(Drug));

    for i = 2:N  % Iterate through Rows  = Z-direciton
        for j = 2:R  % Iterate through Cols = R-direction
            
            % Rb,i* = Net Rate of Drug-CD Binding 
            dDrugCD_dt(i,j) = P1*Drug(i,j)*(P3-Drug_CD(i,j))-Drug_CD(i,j);
            
            % Discretized 2nd order derivative of [Drug] in space 
            % (central difference formula)
            if i > 2 && i < N % If slice is inbetween edge and center 
                d2Drug_dz2(i,j) = (Drug(i+1,j)-2*Drug(i,j)+Drug(i-1,j))/...
                    ((1/N0)^2); % Change in z-direction
            elseif i == 2  % If slice is on center 
                d2Drug_dz2(i,j) = (Drug(i+1,j)-Drug(i,j))/((1/N0)^2);
                % Change in z-direction
            else % i = N % If slice is on edge 
                d2Drug_dz2(i,j) = (-2*Drug(i,j)+Drug(i-1,j))/((1/N0)^2);
                % Change in z-direction
            end 
            
            if j > 2 && j < R % If slice is inbetween edge and center
                d2Drug_dr2(i,j) = (1/j)*(Drug(i,j+1)-Drug(i,j-1))/...
                    (2/R0)+(Drug(i,j+1)-2*Drug(i,j)+Drug(i,j-1))/...
                    ((1/R0)^2); % Change in r-direction
            elseif j == 2  % If slice is on center 
                d2Drug_dr2(i,j) = (1/j)*(Drug(i,j+1)-Drug(i,j))/...
                    (2/R0)+(Drug(i,j+1)-Drug(i,j))/((1/R0)^2); 
                 % Change in r-direction
            else % j == R % If slice is on edge 
                d2Drug_dr2(i,j) = (1/j)*(Drug(i,j)-Drug(i,j-1))/...
                    (2/R0)+(-2*Drug(i,j)+Drug(i,j-1))/((1/R0)^2);
                % Change in r-direction
            end 
            
            % Change in [Drug] in the hydrogel at spatial slice (i) 
                dC_dt(i,j) = P2_z*d2Drug_dz2(i,j) + P2_r*d2Drug_dr2(i,j)...
                    - dDrugCD_dt(i,j);
           
        end
    end 

    % Extract fval as a vector of all other vectors
    % Note: The format of fval should be identical to the input (y)
    % Calculate change in 'Moles Released' in z-dir
    dMR_z_dt = -1/2*P2_z*(Drug(N+1,2:R)-Drug(N,2:R))/(1/N0);
    % Calculate change in 'Moles Released' in r-dir
    dMR_r_dt = -1/2*P2_r*(Drug(2:N,R+1)-Drug(2:N,R))/(1/R0);
    % Output matrix: [ d[Drug]/dt; d[Drug-CD]/dt; d('Moles Released')/dt
    fval = [dC_dt(2:N,2:R); dDrugCD_dt(2:N,2:R); dMR_z_dt; dMR_r_dt'];
    fval = reshape(fval,numel(fval),1); % Restructure output as vector
end
