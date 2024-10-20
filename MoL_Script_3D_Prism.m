% Joshua McGuckin, Samual Huang, Amy Tieu, Julianne Wagner
% Version 03/30/2021

% Function receives input from the GUI for key parameters
function [tSol,cumul_release,output] = MoL_Script_3D_Prism(X,Y,Z,Drug_0,...
    CD_0,K,k2,D,hrs,res,Figure_1,Figure_2,Figure_3,Figure_4,Figure_5,...
    Figure_6,Animation,UploadData)
%% Initialization

% Define hydrogel dimensions
Volume = X*Y*Z; % Hydrogel volume (m^3) 
tSpan = [0 hrs]; % Timespan (hrs)

% Convert timespan into a dimensionless form to reduce computation during
% the ODE solver step 
tau_Span = tSpan * k2; % Dimensionless timespan (t (hrs) * k2 (1/hrs)         
%% Initial Equilibrium concentations

% In order to determine the concentration of free cyclodextrin (Cc) at 
% Equilibrium (t = 0) we solve a binomial for Cc:
% (free CD concentration ie. CD_eq in this code) 
% (free Drug concentration ie. Drug_eq in this code) 
% (Drug-CD complex concentration ie. Drug_CD_eq in this code) 
% General Binomial = a*x^2 + b*x + c = 0  --> 
% eqn = -K*CD_eq^2 - CD_eq - Drug_0*K*CD_eq + CD_0*K*CD_eq + CD_0 = 0; -->
% eqn = -K*(CD_eq^2) + CD_eq*(-Drug_0*K + CD_0*K - 1) + CD_0 = 0; 

% Define the coefficients of the binomial 
a = -K; 
b = -Drug_0*K + CD_0*K - 1 ;
c = CD_0;
equation = [a b c];
r1 = roots(equation);   % Solve the binomial
r_pos = max(r1);        % Extract the positive solution             
                        % (ignore the negative)

% Extract the equilibrium conc. of free Drug, free CD, and Drug-CD complex 
CD_eq = r_pos; % C_C (free CD at the equilibrium (t = 0) )
Drug_CD_eq = Drug_0 - (Drug_0./(1+(K*CD_eq))); % C_LC (equalibrium drug-CD)
Drug_eq = Drug_0./(1+(K*CD_eq)); %Changed from CD_0 (total) to CD_eq (free)                                       
%% Dimensionless Parameters

% Define the dimensionless parameters to simplify computation in the ODE 
% solver step
% P1 = analogous to Drug-to-CD binding affinity 
% P2 = analogous to Diffusivity of drug from CD-conjugated hydrogel (in
% each dimension)
% P3 = analogous to CD-to-Drug conc. ratio (i.e. drug loading)
% For more background info see Fu et al. (2011)

P1 = K*Drug_0;  % K (1/M), Drug_0 (M)
P3 = CD_0/Drug_0; % CD_0 (M), Drug_0 (M)
P2_x = D./(k2*(X./2)^2); % x-dir: D (m^2/hr), k2 (1/hr), h (m)
P2_y = D./(k2*(Y./2)^2); % y-dir: D (m^2/hr), k2 (1/hr), h (m)
P2_z = D./(k2*(Z./2)^2); % z-dir: D (m^2/hr), k2 (1/hr), h (m)
%% Initial Conditions (Put in GUI)

% Resolution of 1D hydrogel space = number of 1D spatial discretizations
X0 = res; % spatial divisions in X-direction
Y0 = res; % spatial divisions in Y-direction
Z0 = res; % spatial divisions in Z-direction

% [Drug](r,t) = [Drug](0,0)
% Center: r = 0 ; Edge: r = R0 
% [Drug](0:edge,0) = Drug_eq --> Drug conc. is uniform at t = 0 
% [Drug-CD](0:edge,0) = Drug_CD_eq
% [CD](0:edge,0) = CD_eq

% Dimensionless Drug and Drug-CD values at each spatial slice
Drug0 = zeros(X0,Y0,Z0); % Preallocate space
Drug0(:) = Drug_eq/Drug_0;  % Uniform [Drug] in all slices 
Drug_CD_0 = zeros(X0,Y0,Z0);
Drug_CD_0(:) = Drug_CD_eq/Drug_0;  % Uniform [Drug-CD] in all slices  
MR_0 = zeros(X0,Y0,Z0);  % No Moles of Drug Released Initally
% Concatenate input arrays [ [Drug]; [Drug-CD]; Moles Released ]
Drug0 = [Drug0 ; Drug_CD_0; MR_0];  % Rows  = Z-direciton; % Cols = R-direction
% Rows = Z-direciton, Cols = R-direction
% Note: The input to the ODE solver needs to be in vector format
y0 = reshape(Drug0,numel(Drug0),1);
%% Solving ODE script

% Plug dimensionless parameters into @domain function and solve ODEs
% Note: ODE15s is recommended, but ODE45 can be used under most conditions
options = odeset('NonNegative',1); % Elimiate negative values (impossible)
[tauSol, ODESol] = ode15s(@(t,y) domain_3D_Prism(t,y,X0,Y0,Z0,P1,P2_x,P2_y,...
    P2_z,P3),tau_Span,y0); 
% DrugSol and tSol are dimensionless

% Convert dimensionless parameters back to original units
ODESol = ODESol*Drug_0; % C_L* = C_L/C_0
tSol = tauSol/k2;           % tau = t * k2
%% Derive Diffusion Contributions in Each Dimension (X,Y,Z)

% Counters 
count = 1;
loc = 1;

% Note: ODESol is a 'numel(tSol) x X0*Y0 x 3' matrix (3D)
for v = 1:(X0*Y0)   % Iteratte through every col of ODESol
    % Extract [Drug] in all directions (x,y,z) and all time points
    DrugSol_vec(:,loc:loc+(X0-1)) = ODESol(:,count:count+(X0-1)); 
    % Extract [Drug-CD] in all directions (x,y,z) and all time points
    DrugCDSol_vec(:,loc:loc+(X0-1)) = ODESol(:,count+X0:count+2*X0-1);
    % Extract 'Moles Released' in both directions (r,z) and all
    % time points
    MRSol_vec(:,loc:loc+(X0-1)) = ODESol(:,count+X0+Y0:count+(X0+Y0+Z0-1)); 
    % Increment counters to next pos.
    count = count + X0+Y0+Z0;
    loc = loc + X0;
end 

% Reshape extracted values back into 3D matrices
DrugSol  = reshape(DrugSol_vec,numel(tSol),X0*Y0,Z0);
DrugCDSol  = reshape(DrugCDSol_vec,numel(tSol),X0*Y0,Z0);
MRSol  = reshape(MRSol_vec,numel(tSol),X0*Y0,Z0);

% Extract [Drug] in hydrogel in X and Y Directions
[~,b,c] = size(DrugSol);
for k = 1:c    % Iterate through all cols in DrugSol matrix (3D)
    for i = 1:numel(tSol)   % Iterate through all timepoints (rows)
        count = 1;  % Counter
        for j = 1:X0 % Iterate through all axial slices (X0)
        % Collect the avg of every X0 rows of DrugSol, shift by X0 
        % each time
            % [Drug] in x-dir
            DrugSol_x(i,j,k) = mean(DrugSol(i,count:(count+X0-1),k));
            % [Drug] in y-dir
            DrugSol_y(i,j,k) = mean(DrugSol(i,j:Y0:b+j-Y0,k));
            % Increment counter to next pos.
            count = count + X0; 
        end 
    end
end

% Extract [Drug] in hydrogel in Z Direction
DrugSol_z = mean(DrugSol,3); % Take the avg of [Drug] matrix in 3rd dim (z)
% Restructure this avg to return avg in the X-Y planes
DrugSol_z = reshape(DrugSol_z,numel(tSol),Z0,Y0); 
% Take avg of the X-Y planes to get [Drug] in z-dir
DrugSol_z = mean(DrugSol_z,2);  

DrugSol_x = mean(DrugSol_x,3); % Find avg of [Drug] in x-dir for XZ planes
DrugSol_y = mean(DrugSol_y,3); % Find avg of [Drug] in y-dir for YZ planes
% Restructure DrugSol_z to same dim as DrugSol_x and DrugSol_y
DrugSol_z = reshape(DrugSol_z,numel(tSol),Z0,1); 
%% Plotting results 

% Input the solved data that represents [Drug] in the hydrogel (DrugSol)
% and other relevant parameters into @Profile to output plots and profile
[cumul_release, output] = profile_3D_Prism(DrugSol,DrugSol_x,DrugSol_y,...
    DrugSol_z,DrugCDSol,tSol,Drug_0,X0,Y0,Z0,X,Y,Z,hrs,Volume,...
    Figure_1,Figure_2,Figure_3,Figure_4,Figure_5,Figure_6,Animation,UploadData);

% Save workspace in current folder
save('Workspace_3D_Prism','CD_0','CD_eq','Drug0','Drug_eq','Drug_0',...
    'Drug_CD_0','Drug_CD_eq','DrugSol','DrugSol_x','DrugSol_y','DrugSol_z',...
    'DrugCDSol','MRSol','X','Y','Z','K','k2','D','hrs','P1','P2_x','P2_y',...
    'P2_z','P3','X0','cumul_release','output','tSol');
end 