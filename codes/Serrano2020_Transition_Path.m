% Optimal Progressivity of Personal Income Tax: A General Equilibrium Evaluation for Spain
% Transition Path between Actual and Optimal Progressivity Levels
% Darío Serrano-Puente
% Banco de España | DG Economics, Statistics and Research | Structural Analysis Division
% darioserrapuente@gmail.com | dario.serrano@bde.es
% https://sites.google.com/view/darioserranopuente/
% August 10th, 2020
%**********************************************************************
% Name: File for acquirement of transition path from baseline to reformed economy
% File name: Serrano2020_Transition_Path.m
% Connected files: in subfolder named "subroutines"

% Description: algorithm to find the transition path between the actual
% economy and an economy with the optimal level of progressivity (and
% average taxes) of the personal income tax.
%**********************************************************************
clc;            % Clear screen
clear;          % Clear memory
rng default;    % Fix the seed for the random number generator (so as to keep simulations the same)
close all;      % Close open windows

addpath(genpath(pwd));
%**********************************************************************

%% 1. Set some options for the code

disp('Running "Optimal Progressivity of Personal Income Tax: A General Equilibrium Evaluation for Spain"')
disp('Transition Path Between Actual and Optimal Economy Algorithm')
disp('Serrano-Puente, Darío (2020)')
disp('Current code version: August 10th, 2020')
disp(' ')
% Parallelize on all possible cores of the computer
delete(gcp('nocreate')); % Close existing interactive parpool session if so
parpool; % Use 2 cores (number of available cores of my computer)
% parpool(1) % Just use 1 core if this line is commented out
PoolDetails=gcp; % Get the number of cores working
Options.NCores=PoolDetails.NumWorkers; % Number of workers/cores to use in the stationary distribution simulation
Options.Tolerance=5*10^(-4); % Iterate until reach this convergence level
Options.Howards=60;
Options.Howards2=500;
Options.Burnin=1000;
Options.Simperiods=10^4; % Simulation periods to get the stationary distribution
Options.NSim=4*10^5; % Used when calculating the statistics on intergenerational earnings correlation and on life-cycle earnings profile. If less than 10^5~10^6, the moments bounce around too much (not enough simulations to really pin it down).
Options.MaxIter=20000;
Options.oldpathweight=0.9; % Set price path to be 'num between 0 and 1' the old path and '1- num between 0 and 1' the new path
Options.weightscheme=4; % Type of update of the price path: ('1') constant weighting, ('2'') exponentially decreasing weighting, ('3') different weighting depending on the 'T', ('4') a mix of cases (2) and (3)
Options.SkipTransition=0; % ('1') to skip transition path finder, ('0') otherwise

disp(' ')
disp('Code options:')
fprintf('Number of cores connected to: %.f \n', Options.NCores)
fprintf('Tolerance level: %.5f \n', Options.Tolerance)
fprintf('Distance of Howards Value Function Iteration Improvement: %.f \n', Options.Howards2)
fprintf('Maximum number of iterations of Howards Value Function Iteration Improvement: %.f \n', Options.Howards)
fprintf('Burnin: %.f \n', Options.Burnin)
fprintf('Simulation periods to get the stationary distribution: %.f \n', Options.Simperiods)
fprintf('Number of simulations for intergenerational correlation and life-cycle profile: %.f \n', Options.NSim)
fprintf('Maximum number of iterations: %.f \n', Options.MaxIter)

%% 2. Load baseline economy (calculated with Calibration algorithm), optimal economy (calculated with Optimal progressivity algorithm), and grids
% These will be the initial and final points on the transition path
 
% 2.1 Load baseline economy and create some objects
load Output/GE/GE_results.mat V Policy StationaryDistr Calib_Params ModelStats_GE ExtraStats_GE
Params=Calib_Params;
ModelStats_BE=ModelStats_GE;
ExtraStats_BE=ExtraStats_GE;
StationaryDistr_BE=StationaryDistr;
V_BE=V;
Policy_BE=Policy;
[e_grid,Gamma,Gamma_ee,~,~]=Create_Exog_Shock(Params); % Gamma transition matrix does not change, then it is sufficient to generate it with baseline parameters

% 2.2 Load baseline economy
load Output/Reform/Opt_results.mat V_opt Policy_opt StationaryDistr_opt Opt_Params ModelStats_opt ExtraStats_opt

% 2.3 Load grids
load Output/Grids.mat d_grid a_grid z_grid n_d n_a n_z
disp(' ')
disp('Grid sizes:')
fprintf('Number of grid points for labor choice: %.f \n', n_d(1))
fprintf('Number of grid points for assets/capital: %.f \n', n_a)
fprintf('Length of state (age) grid: %.f \n', n_z)
N_d=prod(n_d); N_a=prod(n_a); N_z=prod(n_z);

%% 3. Preparation of the transition path finder

% 3.1 Define utility (or return) function equation by calling the return's calculation subroutine
UtilityFn=@(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa) ReturnFn(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa);
UtilityFnParamNames={'r','gamma','varphi','chi','elle','alpha','delta','e1','e2','e3','e4','omega','lambda','tau','kappa'};
DiscountFactorParamNames={'beta'};

% 3.2 Determine next period's assets from current period's decisions
PhiaprimeParamNames={};
Phi_aprimeMatrix=PhiaprimeMatrix(n_d,n_z);

% 3.3 Define GE equations and parameters

% 3.3.1 GE parameters
GEPriceParamNames={'r','omega'}; 
% Note: the government consumption-to-output ratio is fixed at amount 'G',
% then the goverment adjusts its budget by changing the amount devoted to
% transfers to retirees (or people out of the labor market). It means that
% tax revenue gains and losses are redistributed lump-sum to HHs.
p_GE_init=[Params.r,Params.omega]; % Initial values (baseline economy)
p_GE_final=[Opt_Params.r,Opt_Params.omega]; % Final values (optimal economy)

% 3.3.2 Define functions for the General Equilibrium conditions
% The GE equations typically include Market Clearance condtions, but often also things such as Government Budget balance.
% Note: the closer to zero the value given by the function is, the closer the GE condition is to being met.
GeneralEqmEqnParamNames={'G','alpha','delta'}; % Names of parameters that are needed to evaluate the GeneralEqmEqns (parameters not determined as part of GE)

    % (1) Requirement that the interest rate equals the marginal product of capital: i.e. r = alpha*K^(alpha-1)*L^(1-alpha)-delta
GeneralEqmEqn_1 = @(AggVars,p,G,alpha,delta) p(1) - (alpha*(AggVars(1)^(alpha-1))*(AggVars(2)^(1-alpha))-delta); 

    % (2) Government budget balance: i.e. G/Y + Tr/Y = T/Y or G/Y + Tr/Y - T/Y = 0
GeneralEqmEqn_2 = @(AggVars,p,G,alpha,delta) G + AggVars(4)/((AggVars(1)^alpha)*(AggVars(2)^(1-alpha))) - AggVars(3)/((AggVars(1)^alpha)*(AggVars(2)^(1-alpha))); 
% Note: The roles of 'omega', which is contained in p(2), is already captured in the total revenue of income taxes (AggVars(3)) and in total income (Y)).

GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2};

% 3.3.3 Define aggregate SS values
% Used to calculate the integral across the SS distribution
SSvalueParamNames={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
SSvaluesFn_1 = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) a; %K
SSvaluesFn_2 = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) lab*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvaluesFn_TaxRevenue = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa);
SSvaluesFn_Pensions = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) omega*(z>J); % If the agent is retired, he/she earns pension omega (otherwise it is zero).
SSvaluesFn={SSvaluesFn_1,SSvaluesFn_2,SSvaluesFn_TaxRevenue,SSvaluesFn_Pensions};

%% 4. Find the transition path

disp(' ')
disp('Starting calculation of transition path...')

% 4.1 Define T
% This is just a guess and you need to make sure (by checking the output) that the guess is sufficiently large.
T=100;

% 4.2 Define parameter path ('tau', in our case, the change in progressivity parameter)
% If you just want a single change right at the beginning then this is simply the final parameter for all t=1:T. 
% Obviously you can do more advanced stuff like preannoucing changes, or transitory changes.
% ParamPath is of size (T by num_transitionparams)
ParamPath=Opt_Params.tau*ones(T,1);
PathParamNames={'tau'};

% 4.3 Find the equilibrium price path
if Options.SkipTransition == 0
tic
PricePath0=[interp1([1;floor(T/2)],[p_GE_init;p_GE_final],linspace(1,floor(T/2),floor(T/2))'); p_GE_final.*ones(T-floor(T/2),1)];
PricePath=TransitionPath(PricePath0, GEPriceParamNames, ParamPath, PathParamNames, Params, DiscountFactorParamNames, Phi_aprimeMatrix, T, V_opt, StationaryDistr_BE, UtilityFn, UtilityFnParamNames, n_d, n_a, n_z, Gamma, d_grid,a_grid,z_grid, SSvaluesFn,SSvalueParamNames, GeneralEqmEqns, GeneralEqmEqnParamNames,Options.Tolerance,Options.weightscheme,Options.oldpathweight);
runtime3=toc;
save ./SavedOutput/Transition/TransPath.mat PricePath p_GE_init p_GE_final
disp(' ')
fprintf('Finding the equilibrium price path took %.1f seconds. \n', runtime3)
elseif Options.SkipTransition == 1
load ./SavedOutput/Transition/TransPath.mat PricePath p_GE_init p_GE_final
disp(' ')
disp('Equilibrium price path loaded from previous results')
end

% 4.4 Calculate the transition path of other objects
% Now that we have the transition path for prices, let's calculate a bunch of other things associated with it.
% Define the SSvaluesFn_TransPath to give all the aggregate variables we want to see.
%SSvaluesFn_1 = @(d_val,d_ind,a_val,a_ind,s_val,s_ind,pi_s,p_val) a_val; %K (defined above)
%SSvaluesFn_2 = @(d_val,d_ind,a_val,a_ind,s_val,s_ind,pi_s,p_val) d_val(2)*e(s_ind); %L (defined above)
SSvaluesFn_TransPath={SSvaluesFn_1, SSvaluesFn_2};
disp(' ')
disp('ERROR: Aggregate variables path is not implemented yet.')
AggVarsPath=AggVarsPath(SSvaluesFn, SSvalueParamNames,PricePath,PriceParamNames, ParamPath, PathParamNames, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, Phi_aprimeKron_final, Case2_Type,T, V_final, StationaryDist_init, ReturnFn, ReturnFnParamNames, SSvalues_AggVars_final);

save ./SavedOutput/Tranistion/TransPathObjects.mat PricePath AggVarsPath SSvalues_AggVars_init SSvalues_AggVars_final


% Some plots
subplot(2,1,1); plot([p_GE_init;PricePath])
subplot(2,1,2); plot([SSvalues_AggVars_init';AggVarsPath])


