% Joshua McGuckin, Samual Huang, Amy Tieu, Julianne Wagner
% Version 03/30/2021

% Function receives input from the GUI for key parameters
function [tSol,cumul_release,output] = MoL_Script_1D_Slab(r,h,Drug_0,CD_0,...
    K,k2,D,hrs,res,Figure_1,Figure_2,Figure_3,Figure_4,Figure_5,Animation,...
    UploadData)
%% Initialization

% Define hydrogel dimensions
A = r.^2.*pi; % Hydrogel cross-sectional area (m^2)
Volume = A.*h; % Hydrogel volume (m^3) 
tSpan = [0 hrs];  % Timespan (hrs)

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
Drug_eq = Drug_0./(1+(K*CD_eq)); %Changed from CD_0 (total) to CD_eq (free)                                             % concentation [mM]
%% Dimensionless Parameters

% Define the dimensionless parameters to simplify computation in the ODE 
% solver step
% P1 = analogous to Drug-to-CD binding affinity 
% P2 = analogous to Diffusivity of drug from CD-conjugated hydrogel  
% P3 = analogous to CD-to-Drug conc. ratio (i.e. drug loading)
% For more background info see Fu et al. (2011)

P1 = K*Drug_0;  % K (1/M), Drug_0 (M)
P2 = D/(k2*(h/2)^2); % D (m^2/hr), k2 (1/hr), h (m)
P3 = CD_0/Drug_0; % CD_0 (M), Drug_0 (M)
%% Initial Conditions & Spatial Discretizations

% Resolution of 1D hydrogel space = number of 1D spatial discretizations
N0 = res; % res # of 'slices'

% [Drug](z,t) = [Drug](0,0)
% Center: z = 0 ; Edge: z = N0 
% [Drug](0:edge,0) = Drug_eq --> Drug conc. is uniform at t = 0 
% [Drug-CD](0:edge,0) = Drug_CD_eq
% [CD](0:edge,0) = CD_eq

% Dimensionless Drug and Drug-CD values at each spatial slice
% Note: The input to the ODE solver needs to be in vector format
Drug0(1:2:N0*2) = Drug_eq/Drug_0; % Odd rows = [Drug](z,t)
Drug0(2:2:N0*2) = Drug_CD_eq/Drug_0; % Even row = [Drug-CD](z,t)
Drug0(N0*2+1) = 0; % Last row = Moles of Drug released at time(t)
%% Solving ODE script

% Plug dimensionless parameters into @domain function and solve ODEs
% Note: ODE15s is recommended, but ODE45 can be used under most conditions
[tauSol, ODESol] = ode15s(@domain_1D_Slab,tau_Span,Drug0,[],P1,P2,P3); 
% Note: Outputs (tauSol, ODESol) are dimensionless

% Convert dimensionless outputs back to original units
ODESol = ODESol*Drug_0;     % C_L* = C_L/C_0
tSol = tauSol/k2;           % tau = t * k2

% Extract the [Drug](z,t) and [Drug-CD](z,t) from the output ODESol
% This is based on the indexing from the previous section (Drug0 vector)
DrugSol = ODESol(:,1:2:(end-1)); % Odd rows = [Drug](z,t)
DrugCDSol = ODESol(:,2:2:(end-1)); % Even rows = [Drug](z,t)
MRSol = ODESol(:,end); % Last row = Moles of Drug released at time t
%% Plotting results    

% Input the solved data that represents [Drug] in the hydrogel (DrugSol)
% and other relevant parameters into @Profile to output plots and profile
[cumul_release,output] = profile_1D_Slab(DrugSol,tSol,DrugCDSol,Drug_0,...
                    N0,h,r,Volume,Figure_1,Figure_2,Figure_3,...
                    Figure_4,Figure_5,Animation,UploadData);

% Save workspace in current folder
save('Workspace_1D_Slab','CD_0','CD_eq','Drug0','Drug_eq','Drug_0',...
    'Drug_CD_eq','DrugSol','DrugCDSol','MRSol','h','r','K','k2','D','hrs',...
    'P1','P2','P3','N0','cumul_release','output','tSol');
end 
