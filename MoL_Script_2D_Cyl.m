% Joshua McGuckin, Samual Huang, Amy Tieu, Julianne Wagner
% Version 03/30/2021

% Function receives input from the GUI for key parameters
function [tSol,cumul_release,output] = MoL_Script_2D_Cyl(r,h,Drug_0,...
    CD_0,K,k2,D,hrs,res,Figure_1,Figure_2,Figure_3,Figure_4,Figure_5,...
    Figure_6,Animation,UploadData)
%% Initialization

% Define hydrogel dimensions
A = r.^2.*pi; % Hydrogel cross-sectional area (m^2)
Volume = A.*h; % Hydrogel volume (m^3) 
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
P2_z = D/(k2*(h/2)^2); % z-dir: D (m^2/hr), k2 (1/hr), h (m) 
P2_r = D/(k2*(r)^2);  % r-dir: D (m^2/hr), k2 (1/hr), h (m) 
%% Initial Conditions (Put in GUI)

% Resolution of 1D hydrogel space = number of 1D spatial discretizations
N0 = double(res);  % spatial divisions in z-direction
R0 = double(res);  % spatial divisions in r-direction

% [Drug](r,t) = [Drug](0,0)
% Center: r = 0 ; Edge: r = R0 
% [Drug](0:edge,0) = Drug_eq --> Drug conc. is uniform at t = 0 
% [Drug-CD](0:edge,0) = Drug_CD_eq
% [CD](0:edge,0) = CD_eq

% Dimensionless Drug and Drug-CD values at each spatial slice
Drug0 = zeros(N0,R0); % Preallocate space
Drug0(:) = Drug_eq/Drug_0;  % Uniform [Drug] in all slices 
Drug_CD_0 = zeros(N0,R0);   
Drug_CD_0(:) = Drug_CD_eq/Drug_0;  % Uniform [Drug-CD] in all slices 
MR_0 = zeros(2,N0);  % No Moles of Drug Released Initally
% Concatenate input arrays [ [Drug]; [Drug-CD]; Moles Released ]
Drug0 = [Drug0 ; Drug_CD_0; MR_0]; % Rows = Z-direciton, Cols = R-direction
% Note: The input to the ODE solver needs to be in vector format
y0 = reshape(Drug0,numel(Drug0),1); 
%% Solving ODE script

% Plug dimensionless parameters into @domain function and solve ODEs
% Note: ODE15s is recommended, but ODE45 can be used under most conditions
options = odeset('NonNegative',1); % Elimiate negative values (impossible)
[tauSol, ODESol] = ode15s(@domain_2D_Cylinder,tau_Span,y0,[],N0,R0,P1,...
    P2_z,P2_r,P3); 
% DrugSol and tSol are dimensionless

% Convert dimensionless parameters back to original units
ODESol = ODESol*Drug_0; % C_L* = C_L/C_0
tSol = tauSol/k2;           % tau = t * k2
%% Reorganize ODESol into Matrices for [Drug] & [Drug_CD] over Time

% Preallocate space
[n,~] = size(ODESol);
DrugSol = zeros(n,R0*N0);
DrugCDSol = zeros(n,R0*N0);
[a1,b1] = size(DrugSol);
DrugSol_z = zeros(a1,N0);
DrugSol_r = zeros(a1,N0);
MRSol = zeros(a1,2);

for s = 1:a1 % Iteratte through every col and row of ODESol
     numm = 1;  % Counter for raw output from ODEsolver (ODESol)
     numm2 = N0 + R0 + 1;   % Offset for index of 'Moles Released' 
     loc = 1;  % Counter for DrugSol and DrugCDSol matrices
     loc2 = 1;  % Counter for 'Moles Released' matrix
     for z = 1:N0  % Iterate through all 'slices'
         count = numm - 1 + N0; % Increment counter to next pos. of [Drug]
         % Extract [Drug] in both directions (r,z) and all time points
         DrugSol(s,loc:(loc + N0 - 1)) = ODESol(s,numm:count);
         % Extract [Drug-CD] in both directions (r,z) and all time points
         DrugCDSol(s,loc:(loc + N0 - 1)) = ODESol(s,(count + 1):count + N0);
         % Extract 'Moles Released' in both directions (r,z) and all
         % time points
         MRSol(s,loc2:(loc2+1)) = ODESol(s,numm2:numm2+1);
         % Increment counters to next pos.
         numm = numm + N0 + R0 + 2;
         numm2 = numm2 + N0 + R0 + 2;
         loc = loc + N0;
         loc2 = loc2 + 2;
    end 
end 

% Extract [Drug] as a func. of direction (r,z)
for j = 1:a1  % Iterate through all time points (rows)
    for i = 1:N0  % Iterate through all slices (cols)
        % Collect mean of every N0 rows of DrugSol, shift by N0 each time
        % for each axial slice
        count = ((i-1)*N0+1); % Counter
        % Radial [Drug] = mean(row, each block of N0 # idx)value
        DrugSol_z(j,i) = mean(DrugSol(j,count:(count+N0-1)));
        % Radial [Drug] = mean(row, every mth)value
        DrugSol_r(j,i) = mean(DrugSol(j,i:R0:b1-R0+i));
    end
end 
%% Plotting results 

% Input the solved data that represents [Drug] in the hydrogel (DrugSol)
% and other relevant parameters into @Profile to output plots and profile
[cumul_release,output] = profile_2D_Cylinder(DrugSol,DrugSol_r,DrugSol_z,...
    DrugCDSol,tSol,Drug_0,N0,R0,h,r,Volume,Figure_1,Figure_2,Figure_3,...
    Figure_4,Figure_5,Figure_6,Animation,UploadData);

% Save workspace in current folder
save('Workspace_2D_Cylinder','CD_0','CD_eq','Drug0','Drug_eq','Drug_0',...
    'Drug_CD_0','Drug_CD_eq','DrugSol','DrugSol_r','DrugSol_z','DrugCDSol',...
    'MRSol','h','r','K','k2','D','hrs','P1','P2_r','P2_z','P3','N0',...
    'cumul_release','output','tSol');
end 
