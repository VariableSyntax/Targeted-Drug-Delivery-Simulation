% Function is fed with inputs from MoL_Script
% Note: This function is fed with dimensionless paramters!!!

function [fval] = domain_1D_Slab(t,y,P1,P2,P3)
        
    % Define the number of spatial 'slices' 
    N = (numel(y)-1)/2+1; % Step size = 1/N; 
    
    % Define the [Drug] vector
    Drug = zeros(N+1,1); % Preallocate space in a vector
    Drug(2:N) = y(1:2:(length(y)-1)); % Extract all values
                          % [Drug](center + 1:edge - 1,t) from y (input)
    Drug(1) = ((1/3)*(-Drug(3)+4*Drug(2))); %dC_dz (i.e. flux) = 0 at z = 0
    Drug(N+1) = 0;     % [Drug](edge,t) = 0; outside hydrogel is zero
%%
    % Define the [Drug-CD] vector
    Drug_CD = zeros(N+1,1); % Preallocate space in a vector
    Drug_CD(2:N) = y(2:2:(length(y)-1));  % Extract all values
                          % [Drug-CD](center + 1:edge - 1,t) from y (input)
    Drug_CD(1) = ((1/3)*(-Drug_CD(3)+4*Drug_CD(2))); % dCL_dz = 0 at z = 0
    Drug_CD(N+1) = 0;  % [Drug-CD](edge,t) = 0; outside hydrogel is zero
    
    % Define the derivatives
    d2Drug_dz2 = zeros(N+1,1); % Preallocate space in vectors
    dDrugCD_dt = zeros(N+1,1);
    dC_dt = zeros(N+1,1);
    dMR_dt = -1/2*P2*(Drug(N+1)-Drug(N))/(1/N);
    fval = zeros(length(y),1); % Vector to store all other vectors (output)
%%
    for i = 2:N % Iterate throughout all slices between the center and edge 
                 % BUT not the center or edge slices
                 % i = indiv. slice we are looking at
        % Rb,i* = Net Rate of Drug-CD Binding 
        dDrugCD_dt(i) = P1*Drug(i)*(P3-Drug_CD(i))-Drug_CD(i);
        if i > 2 && i < N % If slice i is between the center and edge 
        % Discretized 2nd order derivative (gradient) of [Drug] in 1D space 
        % (central difference formula)
            d2Drug_dz2(i) = (Drug(i+1)-2*Drug(i)+Drug(i-1))/((1/N)^2);
        % Change in [Drug] in the hydrogel at spatial slice (i) 
            dC_dt(i) = P2 * d2Drug_dz2(i) - dDrugCD_dt(i); 
        elseif i == 2 % for slice immediately after the center 
            % (i = 1 is center)
            d2Drug_dz2(i) =(Drug(i+1)-Drug(i))/((1/N)^2);
        % Change in [Drug] in the hydrogel at spatial slice (i) 
            dC_dt(i) = P2*d2Drug_dz2(i)-dDrugCD_dt(i);
        elseif i == N % slice on the edge (i = N in Fu paper)
            d2Drug_dz2(i) =(-2*Drug(i)+Drug(i-1))/((1/N)^2);
        % Change in [Drug] in the hydrogel at spatial slice (i) 
            dC_dt(i) = P2*d2Drug_dz2(i)-dDrugCD_dt(i);
        end
    end
%%
    % Extract fval as a vector of all other vectors
    % Note: The format of fval should be identical to the input (y)
    fval(1:2:length(y)-1) = dC_dt(2:N); % Output d[Drug]/dt
    fval(2:2:length(y)-1) = dDrugCD_dt(2:N); % Output d[Drug-CD]/dt
    fval(end) = dMR_dt; % Output d(Moles of Drug Released)/dt
end
