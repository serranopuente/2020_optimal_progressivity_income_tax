% Optimal Progressivity of Personal Income Tax: A General Equilibrium Evaluation for Spain
% Optimal Progressivity
% Darío Serrano-Puente
% Banco de España | DG Economics, Statistics and Research | Structural Analysis Division
% darioserrapuente@gmail.com | dario.serrano@bde.es
% https://sites.google.com/view/darioserranopuente/
% September 14th, 2020
%**********************************************************************
% Name: File with evaluation of different tax and transfer system reforms
% File name: Serrano2020_Optimal_Progressivity_Reform.m
% Connected files: in subfolder named "subroutines"

% Description: algorithm to find the optimal level of progressivity (and
% average taxes) of the personal income tax. It evaluates a bunch of 
% progressivity reforms and select the optimal one and produce some output.
%**********************************************************************
clc;            % Clear screen
clear;          % Clear memory
rng default;    % Fix the seed for the random number generator (so as to keep simulations the same)
close all;      % Close open windows

addpath(genpath(pwd));
%**********************************************************************

%% 1. Set some options for the code

disp('Running "Optimal Progressivity of Personal Income Tax: A General Equilibrium Evaluation for Spain"')
disp('Optimal Progressivity Finding Algorithm')
disp('Serrano-Puente, Darío (2020)')
disp('Current code version: September 14th, 2020')
disp(' ')
% Parallelize on all possible cores of the computer
delete(gcp('nocreate')); % Close existing interactive parpool session if so
parpool; % Use 2 cores (number of available cores of my computer)
% parpool(1) % Just use 1 core if this line is commented out
PoolDetails=gcp; % Get the number of cores working
Options.NCores=PoolDetails.NumWorkers; % Number of workers/cores to use in the stationary distribution simulation
Options.Guess=1; % Make initial guess for optimal progressivity finding loop from previous try ('1') or no guess ('0'). Note: make sure that the dimensions of the grid of the previous try match the dimensions of the selected grid for evaluation
Options.SkipLoop=1; % Skip progressivity loop ('1') or activate it ('0')
Options.SkipReformFind=0; % Skip part of progressivity loop devoted to find the equilibrium of the reformed economy ('1') or activate it ('0')
Options.SkipLoop_decomp=1; % Skip loop to decompose welfare changes ('1') or activate it ('0')
Options.SkipLoop_HHtype=1; % Skip loop to get welfare changes by household type ('1') or activate it ('0')
Options.Guess_HHtype=1; % Make initial guess for welfare changes by household type from previous try ('1') or no guess ('0'). Note: make sure that the dimensions of the grid of the previous try match the dimensions of the selected grid for evaluation
Options.Tolerance=10^(-5); % Iterate until reach this convergence level
Options.fminoptions = optimset('TolFun',10^(-3),'FunValCheck','on','MaxIter',150); % Options for minimum search algorithm - The tolerance is downgraded a bit, because we are checking a quick overview of where the optimal progressivity value could lie in terms of welfare gains (CEV)
Options.Howards=60;
Options.Howards2=500;
Options.Burnin=1000;
Options.Simperiods=10^5; % Simulation periods to get the stationary distribution
Options.NSim=8*10^5; % Used when calculating the statistics on intergenerational earnings correlation and on life-cycle earnings profile. If less than 10^5~10^6, the moments bounce around too much (not enough simulations to really pin it down).
Options.MaxIter=40000;

disp(' ')
disp('Code options:')
fprintf('Number of cores connected to: %.f \n', Options.NCores)
fprintf('Retriving guess for optimal progressivity finder from previous tries: %.f \n', Options.Guess)
fprintf('Skik big loop of progressivity reforms: %.f \n', Options.SkipLoop)
fprintf('Tolerance level: %.5f \n', Options.Tolerance)
fprintf('Distance of Howards Value Function Iteration Improvement: %.f \n', Options.Howards2)
fprintf('Maximum number of iterations of Howards Value Function Iteration Improvement: %.f \n', Options.Howards)
fprintf('Burnin: %.f \n', Options.Burnin)
fprintf('Simulation periods to get the stationary distribution: %.f \n', Options.Simperiods)
fprintf('Number of simulations for intergenerational correlation and life-cycle profile: %.f \n', Options.NSim)
fprintf('Maximum number of iterations: %.f \n', Options.MaxIter)

%% 2. Load baseline economy and grids (calculated with Calibration algorithm)
 
% 2.1 Load baseline economy and create some objects
load Output/GE/GE_results.mat p_eqm MarketClearance V Policy StationaryDistr Calib_Params ModelStats_GE ExtraStats_GE GEtime
Params=Calib_Params;
ModelStats_BE=ModelStats_GE;
ExtraStats_BE=ExtraStats_GE;
StationaryDistr_BE=StationaryDistr;
V_BE=V;
Policy_BE=Policy;
[e_grid,Gamma,Gamma_ee,gammastar,gammastarfull]=Create_Exog_Shock(Params);

% 2.2 Load grids
load Output/Grids.mat d_grid a_grid z_grid n_d n_a n_z
disp(' ')
disp('Grid sizes:')
fprintf('Number of grid points for labor choice: %.f \n', n_d(1))
fprintf('Number of grid points for assets/capital: %.f \n', n_a)
fprintf('Length of state (age) grid: %.f \n', n_z)
N_d=prod(n_d); N_a=prod(n_a); N_z=prod(n_z);

%% 3. Preparation of progressivity reform loop

% 3.1 Create grid for possible values of progressivity to find the optimal one
tau_grid=[0.00:0.010:0.30,0.32:0.02:0.34,0.37:0.03:0.40,0.45:0.05:0.50]; % Sparse grid concentrated at maximum
MaxReforms=length(tau_grid);
fprintf('Number of grid points for progressivity level to find the optimal one: %.f \n', MaxReforms)

% 3.2 Define utility (or return) function equation by calling the return's calculation subroutine
UtilityFn=@(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa) ReturnFn(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa);
UtilityFn_CEV=@(CEV,lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa) ReturnFn_CEV(CEV,lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa); % Same as UtilityFn but flexible to evaluate CEV to compute welfare changes in % of consumption
UtilityFnParamNames={'r','gamma','varphi','chi','elle','alpha','delta','e1','e2','e3','e4','omega','lambda','tau','kappa'};
DiscountFactorParamNames={'beta'};

% 3.3 Determine next period's assets from current period's decisions
PhiaprimeParamNames={};
Phi_aprimeMatrix=PhiaprimeMatrix(n_d,n_z);

% 3.4 Create counter of the possible reforms of the progressivity
ReformCounter=0;
save ./Output/Reform/ReformCounter.mat ReformCounter MaxReforms

% 3.5 Define parameters that would change to define the new equilibrium and
GEPriceParamNames={'r','lambda','omega'}; 
% Note: the government consumption-to-output ratio is fixed at amount 'G/Y',
% then the goverment adjusts its budget by changing the level of 'lambda' 
% as response to change in progressivty ('tau') and 'omega' to mantain Tr/Y ratio.
TargetsBudget=ModelStats_BE(2:3); % Level of G/Y and Tr/Y of the baseline economy

% 3.6 Augment weight to government budget balance condition in GE of the
% reformed economies, otherwise this target would be never achieved
weight_BudgetBalance = [10 20]';

% 3.7 Guesses and initial values
V0=reshape(ones(n_a,n_z),[N_a,N_z]);
Ref_Params=Params; % In order not to overwrite parameters of the baseline economy
if Options.Guess == 1
% Take the initial values for 'r' and 'lambda' (parameters that solve the equilibrium of the reformed economy) arbitrarily or from previous tries
% Note: make sure that the dimensions of the grid of the previous try match the dimensions of the selected grid for evaluation
load ./Output/Reform/OutputRef.mat OutputRef CovergenceAlgorithms
OutputRef_Old = OutputRef;
CovergenceAlgorithms_Old = CovergenceAlgorithms;
r_init = 0.0820;
lambda_init = 1.36; 
omega_init = 2.53;
% Give an initial guess on CEV for the aggregate welfare change computation arbitrarily or from previous tries
CEV_init=OutputRef_Old(39,1:size(OutputRef_Old,2))/100;
end

% 3.8 Generate matrices to store the results of the loop
OutputRef=nan(39,MaxReforms);
CovergenceAlgorithms=nan(2,MaxReforms);

%%  4. Progressivity reform loop
if Options.SkipLoop == 0
disp(' ')
disp('Initializing progressivity reform loop...')

for jj = 1:MaxReforms
if Options.SkipReformFind == 0
% 4.1 Find the equilibrium of the reformed economy
load ./Output/Reform/ReformCounter.mat ReformCounter MaxReforms
ReformCounter=ReformCounter+1;
save ./Output/Reform/ReformCounter.mat ReformCounter MaxReforms
IterReformCounter=0;
save ./Output/Reform/IterReformCounter.mat IterReformCounter
if Options.Guess == 1
% Start minimum search with previous tries findings
Ref_Params.r=r_init;
Ref_Params.lambda=lambda_init;
Ref_Params.omega=omega_init;
CEV_init_scalar=CEV_init(jj);
end
disp(' ')
fprintf('Intializing General Equilibrium calculation of reform %.f out of %.f... \n', jj, MaxReforms);
tic;
ReformEqmFn = @(p) ReformFn(p, tau_grid(jj), weight_BudgetBalance, Ref_Params, DiscountFactorParamNames, UtilityFn, UtilityFnParamNames, GEPriceParamNames, V0, Phi_aprimeMatrix, TargetsBudget, PhiaprimeParamNames, n_d, n_a, n_z, a_grid, d_grid, z_grid, Gamma, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter);
p0=nan(length(GEPriceParamNames),1); % Identify parameters that are solvers of the equilibrium to get an initial guess (which is the value that these parameters take in the baseline economy)
for ii=1:length(GEPriceParamNames)
    p0(ii)=Ref_Params.(GEPriceParamNames{ii});
end
[p_eqm_ref,fval]=fminsearch(ReformEqmFn,p0,Options.fminoptions); % Search for the values of the parameters that solve the equilibrium
CovergenceAlgorithms(1,jj)=fval;
GEref_time=toc;
disp(' ')
fprintf('End of General Equilibrium calculation of reform %.f out of %.f. Time used was: %.f. \n', jj, MaxReforms, GEref_time);

% Load previous reformed economy equilibrium if already found in previous runs of the loop
elseif Options.SkipReformFind == 1
disp(' ')
fprintf('General Equilibrium from reformed economy retrieved from previous runs of the loop. \n')
Ref_Params.r=OutputRef_Old(2,jj);
Ref_Params.lambda=OutputRef_Old(3,jj);
Ref_Params.omega=OutputRef_Old(4,jj);
CovergenceAlgorithms(1,jj)=CovergenceAlgorithms_Old(1,jj);
CEV_init_scalar=CEV_init(jj);
end

% 4.2 Get VF and stationary distribution of the new equilibrium of the reformed economy
disp(' ')
fprintf('Evaluating a few objects at General Equilibrium (getting VF and stationary distribution of the reformed economy %.f out of %.f)... \n', jj, MaxReforms);
if Options.SkipReformFind == 0
Ref_Params.r=p_eqm_ref(1);
r_init=p_eqm_ref(1);
Ref_Params.lambda=p_eqm_ref(2);
lambda_init=p_eqm_ref(2);
Ref_Params.omega=p_eqm_ref(3);
omega_init=p_eqm_ref(3);
end
Ref_Params.tau=tau_grid(jj);
Ref_Params.w=(1-Ref_Params.alpha)*(((Ref_Params.r+Ref_Params.delta)/(Ref_Params.alpha))^(Ref_Params.alpha/(Ref_Params.alpha-1)));
[V_ref, Policy_ref]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Ref_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2);
StationaryDistr_ref=StationaryDist(Policy_ref,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,Options.NCores,Options.Simperiods,Options.Burnin,Options.Tolerance,Options.MaxIter);
% Save the output of the parameters of the reform
OutputRef(1,jj)=Ref_Params.tau; % Progressivity value ('tau')
OutputRef(2,jj)=Ref_Params.r; % Interest rate of the new equilibrium following the reform ('r')
OutputRef(3,jj)=Ref_Params.lambda; % Average level of taxes in HSV tax funcion specification ('lambda')
OutputRef(4,jj)=Ref_Params.omega; % Normalized transfers to retirees ('omega')
OutputRef(5,jj)=Ref_Params.w; % Wage level after reform in new equilibrium ('w')

% 4.3 Compute some measures/statistics and save the output
disp(' ')
fprintf('Retrieving model statistics of the reformed economy %.f out of %.f... \n', jj, MaxReforms);
SSvalueParamNames=struct();
SSvalueParamNames(1).Names={};
SSvaluesFn_K = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
SSvaluesFn_L = @(lab,anext,a,z,e1,e2,e3,e4) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvalueParamNames(3).Names={};
SSvaluesFn_H = @(lab,anext,a,z) lab; %H
SSvalueParamNames(4).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
SSvaluesFn_TaxRevenue = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa);
SSvalueParamNames(5).Names={'J','omega'};
SSvaluesFn_Pensions = @(lab,anext,a,z,J,omega) omega*(z>J); % If the agent is retired, he/she earns pension omega (otherwise it is zero).
SSvalueParamNames(6).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'}; % SSvaluesFn_Consumption
SSvaluesFn_Consumption = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) ConsumptionFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa); % C
SSvalueParamNames(7).Names={'delta'};
SSvaluesFn_I = @(lab,anext,a,z,delta) anext-a*(1-delta); %I
SSvaluesFn={SSvaluesFn_K,SSvaluesFn_L,SSvaluesFn_H,SSvaluesFn_TaxRevenue,SSvaluesFn_Pensions,SSvaluesFn_Consumption,SSvaluesFn_I};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDistr_ref, Policy_ref, SSvaluesFn, Ref_Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);

% 4.3.1 Macroeconomic aggregates and factor ratios
OutputRef(6,jj)=(SSvalues_AggVars(1)^Ref_Params.alpha)*(SSvalues_AggVars(2)^(1-Ref_Params.alpha)); % Y (Production technology function)
OutputRef(7,jj)=SSvalues_AggVars(1); % K
OutputRef(8,jj)=SSvalues_AggVars(2); % L 
OutputRef(9,jj)=SSvalues_AggVars(6); % C
OutputRef(10,jj)=(SSvalues_AggVars(3)/Ref_Params.elle)*100; % H/elle (share of disposable time allocated to market activities)
OutputRef(11,jj)=(SSvalues_AggVars(1)/SSvalues_AggVars(2)); % K/L
OutputRef(12,jj)=(SSvalues_AggVars(2)/SSvalues_AggVars(3)); % L/H
OutputRef(13,jj)=(OutputRef(6,jj)/SSvalues_AggVars(3)); % Y/H
OutputRef(14,jj)=(SSvalues_AggVars(6)/SSvalues_AggVars(3)); % C/H
OutputRef(15,jj)=(SSvalues_AggVars(1)/OutputRef(6,jj)); % K/Y

% 4.3.2 Expenditure and fiscal policy ratios
OutputRef(16,jj)=(SSvalues_AggVars(6)/OutputRef(6,jj))*100; % C/Y
OutputRef(17,jj)=(SSvalues_AggVars(7)/OutputRef(6,jj))*100; % I/Y
OutputRef(18,jj)=((SSvalues_AggVars(4)-SSvalues_AggVars(5))/OutputRef(6,jj))*100; % G/Y (Government-expenditure to output ratio, G=T-Tr=Tax Revenue - Pensions)
OutputRef(19,jj)=(SSvalues_AggVars(5)/OutputRef(6,jj))*100; % Tr/Y
OutputRef(20,jj)=(SSvalues_AggVars(4)/OutputRef(6,jj))*100; % T/Y

% 4.3.3 Ratio of earnings of old-to-young and intergenerational correlation
e=[Ref_Params.e1,Ref_Params.e2,Ref_Params.e3,Ref_Params.e4,0,0,0,0];
OutputRef(21,jj)=RatioEarningsOldYoung(Options.NSim, StationaryDistr_ref, Policy_ref, Phi_aprimeMatrix, n_d,n_a,n_z, d_grid, Gamma, e,Ref_Params.w,Ref_Params.J); % Ratio of earnings of 40 year olds to 20 year olds
OutputRef(22,jj)=IntergenerationalEarnings(Options.NSim,StationaryDistr_ref, Policy_ref, Phi_aprimeMatrix, n_d,n_a,n_z,d_grid, Gamma, e,Ref_Params.w,Ref_Params.J); % Cross sectional correlation of incomes between fathers and sons

% 4.3.4 Distributional statistics of income and wealth
SSvalueParamNames(1).Names={'e1','e2','e3','e4','w','r','J','omega'}; % income
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J);
SSvalueParamNames(2).Names={}; % K
SSvaluesFn={SSvaluesFn_income,SSvaluesFn_K};
SSvalues_LorenzCurves=SSvalues_LorenzCurve(StationaryDistr_ref, Policy_ref, SSvaluesFn, Ref_Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,100);
OutputRef(23,jj)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(1,:)); % Gini for Income
OutputRef(24:27,jj)=100*(SSvalues_LorenzCurves(1,[40,60,80,100])-SSvalues_LorenzCurves(1,[1,41,61,81])); % Income Lorenz Curve - Quintiles (%) 
OutputRef(28:30,jj)=100*(SSvalues_LorenzCurves(1,[95,99,100])-SSvalues_LorenzCurves(1,[90,95,99])); % Income Lorenz Curve - 90-95, 95-99, and 99-100 (%)
OutputRef(31,jj)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(2,:)); % Gini for Wealth
OutputRef(32:35,jj)=100*(SSvalues_LorenzCurves(2,[40,60,80,100])-SSvalues_LorenzCurves(2,[1,41,61,81])); % Wealth Lorenz Curve - Quintiles (%)
OutputRef(36:38,jj)=100*(SSvalues_LorenzCurves(2,[95,99,100])-SSvalues_LorenzCurves(2,[90,95,99])); % Wealth Lorenz Curve - 90-95, 95-99, and 99-100 (%)
% Adjust income and wealth shares to make them to sum up to 1 - this is due
% to calculus error that mak it impossible to sum up to 1 directly (it
% normally results in a sum of 98 instead of 100)
TotIncomeShares=sum(OutputRef(24:27,jj));
OutputRef(24:27,jj)=(OutputRef(24:27,jj))/TotIncomeShares*100;
TotWealthShares=sum(OutputRef(32:35,jj));
OutputRef(32:35,jj)=(OutputRef(32:35,jj))/TotWealthShares*100;

% 4.3.5 Compute aggregate welfare changes (CEV) derived from the reform
IterCEVCounter=0;
save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
AggFVDist = @(CEV) CompareVF(V0, CEV, V_ref, StationaryDistr_ref, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn_CEV, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter, 0);
[AggWelfareChange,fval]=fminsearch(AggFVDist,CEV_init_scalar,Options.fminoptions); % CEV initial resulting from previous try
CovergenceAlgorithms(2,jj)=fval;
% CEV_init=AggWelfareChange; % Update guess on CEV for next reformed economy with the Aggregate Welfare change calculated in this reformed economy
OutputRef(39,jj)=AggWelfareChange*100;
disp(' ')
fprintf('End of calculating CEV to determine welfare gains of the reformed economy %.f out of %.f. \n', jj, MaxReforms);

% 4.3.6 Save results on disk
save ./Output/Reform/OutputRef.mat OutputRef CovergenceAlgorithms

disp(' ')
fprintf('End of calculation of GE and statistics of the reform %.f out of %.f. \n', jj, MaxReforms);
end
disp(' ')
disp('End of progressivity reform loop.')

% 4.4 Load previously calculated results if progressivity reform loop is not activated
elseif Options.SkipLoop == 1
    disp(' ')
    disp('Progressivity reform loop is not activated. Previous results are loaded.')
    load ./Output/Reform/OutputRef.mat OutputRef CovergenceAlgorithms
end

% 4.5 Include statistics / measures of the baseile economy
 disp(' ')
 disp('Computing statistics of the baseline economy, such that they are comparable with the reformed economy ones...')
% 4.5.1 Save parameters
OutputRef_actual=nan(39,MaxReforms+1); OutputRef_actual(:,1:MaxReforms)=OutputRef;
OutputRef_actual(1,MaxReforms+1)=Params.tau; % Progressivity value ('tau')
OutputRef_actual(2,MaxReforms+1)=Params.r; % Interest rate ('r')
OutputRef_actual(3,MaxReforms+1)=Params.lambda; % Average level of taxes in HSV tax funcion specification ('lambda')
OutputRef_actual(4,MaxReforms+1)=Params.omega; % Normalized transfers to retirees ('omega')
OutputRef_actual(5,MaxReforms+1)=Params.w; % Wage level ('w')

% 4.5.2 Compute some measures/statistics and save the output
SSvalueParamNames=struct();
SSvalueParamNames(1).Names={};
SSvaluesFn_K = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
SSvaluesFn_L = @(lab,anext,a,z,e1,e2,e3,e4) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvalueParamNames(3).Names={};
SSvaluesFn_H = @(lab,anext,a,z) lab; %H
SSvalueParamNames(4).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
SSvaluesFn_TaxRevenue = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa);
SSvalueParamNames(5).Names={'J','omega'};
SSvaluesFn_Pensions = @(lab,anext,a,z,J,omega) omega*(z>J); % If the agent is retired, he/she earns pension omega (otherwise it is zero).
SSvalueParamNames(6).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'}; % SSvaluesFn_Consumption
SSvaluesFn_Consumption = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) ConsumptionFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa); % C
SSvalueParamNames(7).Names={'delta'};
SSvaluesFn_I = @(lab,anext,a,z,delta) anext-a*(1-delta); %I
SSvaluesFn={SSvaluesFn_K,SSvaluesFn_L,SSvaluesFn_H,SSvaluesFn_TaxRevenue,SSvaluesFn_Pensions,SSvaluesFn_Consumption,SSvaluesFn_I};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDistr_BE, Policy_BE, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);

% 4.5.3 Macroeconomic aggregates and factor ratios
OutputRef_actual(6,MaxReforms+1)=(SSvalues_AggVars(1)^Params.alpha)*(SSvalues_AggVars(2)^(1-Params.alpha)); % Y (Production technology function)
OutputRef_actual(7,MaxReforms+1)=SSvalues_AggVars(1); % K
OutputRef_actual(8,MaxReforms+1)=SSvalues_AggVars(2); % L 
OutputRef_actual(9,MaxReforms+1)=SSvalues_AggVars(6); % C
OutputRef_actual(10,MaxReforms+1)=(SSvalues_AggVars(3)/Params.elle)*100; % H/elle (share of disposable time allocated to market activities)
OutputRef_actual(11,MaxReforms+1)=(SSvalues_AggVars(1)/SSvalues_AggVars(2)); % K/L
OutputRef_actual(12,MaxReforms+1)=(SSvalues_AggVars(2)/SSvalues_AggVars(3)); % L/H
OutputRef_actual(13,MaxReforms+1)=(OutputRef_actual(6,MaxReforms+1)/SSvalues_AggVars(3)); % Y/H
OutputRef_actual(14,MaxReforms+1)=(SSvalues_AggVars(6)/SSvalues_AggVars(3)); % C/H
OutputRef_actual(15,MaxReforms+1)=(SSvalues_AggVars(1)/OutputRef_actual(6,MaxReforms+1)); % K/Y

% 4.5.4 Expenditure and fiscal policy ratios
OutputRef_actual(16,MaxReforms+1)=(SSvalues_AggVars(6)/OutputRef_actual(6,MaxReforms+1))*100; % C/Y
OutputRef_actual(17,MaxReforms+1)=(SSvalues_AggVars(7)/OutputRef_actual(6,MaxReforms+1))*100; % I/Y
OutputRef_actual(18,MaxReforms+1)=((SSvalues_AggVars(4)-SSvalues_AggVars(5))/OutputRef_actual(6,MaxReforms+1))*100; % G/Y (Government-expenditure to output ratio, G=T-Tr=Tax Revenue - Pensions)
OutputRef_actual(19,MaxReforms+1)=(SSvalues_AggVars(5)/OutputRef_actual(6,MaxReforms+1))*100; % Tr/ Y
OutputRef_actual(20,MaxReforms+1)=(SSvalues_AggVars(4)/OutputRef_actual(6,MaxReforms+1))*100; % T/Y

% 4.5.5 Ratio of earnings of old-to-young and intergenerational correlation
e=[Params.e1,Params.e2,Params.e3,Params.e4,0,0,0,0];
OutputRef_actual(21,MaxReforms+1)=RatioEarningsOldYoung(Options.NSim, StationaryDistr_BE, Policy_BE, Phi_aprimeMatrix, n_d,n_a,n_z, d_grid, Gamma, e,Params.w,Params.J); % Ratio of earnings of 40 year olds to 20 year olds
OutputRef_actual(22,MaxReforms+1)=IntergenerationalEarnings(Options.NSim,StationaryDistr_BE, Policy_BE, Phi_aprimeMatrix, n_d,n_a,n_z,d_grid, Gamma, e,Params.w,Params.J); % Cross sectional correlation of incomes between fathers and sons

% 4.5.6 Distributional statistics of income and wealth
SSvalueParamNames(1).Names={'e1','e2','e3','e4','w','r','J','omega'}; % income
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J);
SSvalueParamNames(2).Names={}; % K
SSvaluesFn={SSvaluesFn_income,SSvaluesFn_K};
SSvalues_LorenzCurves=SSvalues_LorenzCurve(StationaryDistr_BE, Policy_BE, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,100);
OutputRef_actual(23,MaxReforms+1)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(1,:)); % Gini for Income
OutputRef_actual(24:27,MaxReforms+1)=100*(SSvalues_LorenzCurves(1,[40,60,80,100])-SSvalues_LorenzCurves(1,[1,41,61,81])); % Income Lorenz Curve - Quintiles (%) 
OutputRef_actual(28:30,MaxReforms+1)=100*(SSvalues_LorenzCurves(1,[95,99,100])-SSvalues_LorenzCurves(1,[90,95,99])); % Income Lorenz Curve - 90-95, 95-99, and 99-100 (%)
OutputRef_actual(31,MaxReforms+1)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(2,:)); % Gini for Wealth
OutputRef_actual(32:35,MaxReforms+1)=100*(SSvalues_LorenzCurves(2,[40,60,80,100])-SSvalues_LorenzCurves(2,[1,41,61,81])); % Wealth Lorenz Curve - Quintiles (%)
OutputRef_actual(36:38,MaxReforms+1)=100*(SSvalues_LorenzCurves(2,[95,99,100])-SSvalues_LorenzCurves(2,[90,95,99])); % Wealth Lorenz Curve - 90-95, 95-99, and 99-100 (%)
% Adjust income and wealth shares to make them to sum up to 1 - this is due
% to calculus error that mak it impossible to sum up to 1 directly (it
% normally results in a sum of 98 instead of 100)
TotIncomeShares=sum(OutputRef_actual(24:27,MaxReforms+1));
OutputRef_actual(24:27,MaxReforms+1)=(OutputRef_actual(24:27,MaxReforms+1))/TotIncomeShares*100;
TotWealthShares=sum(OutputRef_actual(32:35,MaxReforms+1));
OutputRef_actual(32:35,MaxReforms+1)=(OutputRef_actual(32:35,MaxReforms+1))/TotWealthShares*100;

% 4.5.7 Compute aggregate welfare changes (CEV) are 0, since it is the baseline economy
OutputRef_actual(39,MaxReforms+1)=0;

% 4.6 Save results on disk and export data for graphs in Excel
save ./Output/Reform/OutputRef_actual.mat OutputRef_actual
xlswrite('Output\Tables\OutputRef_actual.xlsx',OutputRef_actual)

%% 5. Find optimal progressivity (and transfers) reform

% Select the combination of progressivity ('tau') and average level of
% taxes ('lambda') that maximizes the aggregate welfare changes (CEV)
% keeping G/Y and Tr/Y constant and having Gov. Budget Balance

% 5.1 Sort and get parameters of the optimal reform
OutputRef_actual=sortrows(OutputRef_actual',39)'; % Transpose matrix to be able to sort it entirely (sort by Aggregate Welfare change)
Opt_Params=Params;
Opt_Params.tau=OutputRef_actual(1,size(OutputRef_actual,2));   % Progressivity value ('tau')
Opt_Params.r=OutputRef_actual(2,size(OutputRef_actual,2));     % Interest rate of the new equilibrium following the optimal reform ('r')
Opt_Params.lambda=OutputRef_actual(3,size(OutputRef_actual,2));% Average level of taxes in HSV tax funcion specification ('lambda')
Opt_Params.omega=OutputRef_actual(4,size(OutputRef_actual,2)); % Normalized transfers to retirees implied by the optimal reform ('omega')
Opt_Params.w=OutputRef_actual(5,size(OutputRef_actual,2));     % Wages implied by optimal reform ('w')

% 5.2 Get VF and stationary distribution of the optimal reform
[V_opt, Policy_opt]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Opt_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2);
StationaryDistr_opt=StationaryDist(Policy_opt,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,Options.NCores,Options.Simperiods,Options.Burnin,Options.Tolerance,Options.MaxIter);

% 5.3 Get model statistics of optimal reform
ModelStats_opt=ModelStatistics(Options.NSim,StationaryDistr_opt, Policy_opt, Phi_aprimeMatrix, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Opt_Params);
ExtraStats_opt=ModelStatistics_ExtraStats(StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params);

% 5.4 Save optimal reform results
save Output/Reform/Opt_results.mat V_opt Policy_opt StationaryDistr_opt Opt_Params ModelStats_opt ExtraStats_opt

%% 6. Decompose aggregate welfare changes (CEV) derived from the optimal reform

if Options.SkipLoop_decomp == 0
% Generate two ficticious economies to be able to decomposte the aggregate welfare changes
disp(' ')
disp('Starting loop to decompose aggregate welfare changes...')
AggWelfareChange_decomp=nan(4,1);
AggWelfareChange_decomp(1,1) = OutputRef_actual(39,size(OutputRef_actual,2)); % Aggregate welfare change
CovergenceAlgorithms_decomp=nan(2,1);
CEV_init=0; % Initial guess is set to 0.
ReformCounter=1;
MaxReforms=1;

% 6.1 Economy ignoring changes in the distribution of HH (use stationary distribution of the baseline economy)
IterCEVCounter=0;
save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
AggFVDist_a = @(CEV) CompareVF(V0, CEV, V_opt, StationaryDistr_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn_CEV, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter, 1);
[AggWelfareChange_a,fval]=fminsearch(AggFVDist_a,CEV_init,Options.fminoptions);
CovergenceAlgorithms_decomposition(1,1)=fval;
disp(' ')
fprintf('End of calculating CEV to determine welfare gains of a hypothetical economy ignoring changes in the distribution of HH after reform %.f out of %.f. \n', ReformCounter, MaxReforms);
 
% 6.2 Economy ignoring changes in the distribution of HH and changes in size of the economy/equilibrium prices (use prices and stationary distribution of the baseline economy, the only different thing is the tax & transfer system  with respect to the baseline economy)
IterCEVCounter=0;
save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
b_Params=Params;
b_Params.tau=Opt_Params.tau;
b_Params.lambda=Opt_Params.lambda;
b_Params.omega=Opt_Params.omega;
[V_b, Policy_b]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, b_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2);
AggFVDist_b = @(CEV) CompareVF(V0, CEV, V_b, StationaryDistr_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn_CEV, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter, 2);
[AggWelfareChange_b,fval]=fminsearch(AggFVDist_b,CEV_init,Options.fminoptions);
CovergenceAlgorithms_decomposition(2,1)=fval;
disp(' ')
fprintf('End of calculating CEV to determine welfare gains of a hypothetical economy ignoring changes in the distribution of HH and changes in size of the economy/equilibrium prices after reform %.f out of %.f. \n', ReformCounter, MaxReforms);
 
% 6.3 Calculate and save share of the aggregate welfare change (CEV) due to each motive and aggregate welfare change (CEV) itself
AggWelfareChange_TT = (AggWelfareChange_b*100/AggWelfareChange_decomp(1,1)); % Aggregate welfare changes due to reshiffling resources between HHs (change in tax & transfers system)
AggWelfareChange_GE = ((AggWelfareChange_a*100 - AggWelfareChange_b*100)/AggWelfareChange_decomp(1,1)); % Aggregate welfare changes due to GE effects only (change in equilibrium prices / size of the economy)
AggWelfareChange_Dist = ((AggWelfareChange_decomp(1,1) - AggWelfareChange_a*100)/AggWelfareChange_decomp(1,1)); % Aggregate welfare changes due to change in distribution of HHs (change in stationary distribution)
AggWelfareChange_decomp(2,1) = AggWelfareChange_TT;
AggWelfareChange_decomp(3,1) = AggWelfareChange_GE;
AggWelfareChange_decomp(4,1) = AggWelfareChange_Dist;
disp(' ')
disp('End of loop to decompose aggregate welfare changes.')

% 6.4 Save results on disk or load previous ones
save Output/Reform/AggWelfareChange_decomp.mat AggWelfareChange_decomp

elseif Options.SkipLoop_decomp == 1
    disp(' ')
    disp('Loop to decompose aggregate welfare changes is not activated. Previous results are loaded.')
    load Output/Reform/AggWelfareChange_decomp.mat AggWelfareChange_decomp
end

%% 7. Compute welfare changes (CEV) by type of household in the optimal reform scenario

if Options.SkipLoop_HHtype == 0
disp(' ')
disp('Starting loop to compute welfare changes (CEV) by type of household (ranked by income, wealth, and VF).')
% Note: we want to classify the households by wealth, by income, and by value function
CEV_init=0; % Initial guess is set to 0.

% Uncomment following lines to get decomposition by HH type sorted by income & wealth
% 7.1 Sort by wealth
n_pctile_groups = 10;
pctile.lb = (0:(100/n_pctile_groups):100-(100/n_pctile_groups))';
pctile.ub = pctile.lb+(100/n_pctile_groups);
pctile.ub(length(pctile.ub))=pctile.ub(length(pctile.ub))+0.00001; % To capture the very right tail of the distribution

if Options.Guess_HHtype == 1
    load Output/Reform/WelfareChange_HH_wealth.mat WelfareChange_HH_wealth CovergenceAlgorithms_HH_wealth
elseif Options.Guess_HHtype == 0
    WelfareChange_HH_wealth = ones(11,n_pctile_groups);
    CovergenceAlgorithms_HH_wealth = ones(9,n_pctile_groups);
end

SSvalueParamNames=struct();
SSvalueParamNames(1).Names={};
SSvaluesFn_K = @(lab,anext,a,z) a; %K
SSvaluesFn={SSvaluesFn_K};

for ww = 1:length(pctile.lb)
    for ss = 0:8 % 0 denotes all states
        if CovergenceAlgorithms_HH_wealth(1+ss,ww) > 10^(-3)
            IterCEVCounter=0;
            save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
            disp(' ')
            fprintf('Calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by wealth for income shock s = %.f... \n', pctile.lb(ww), pctile.ub(ww), ss);
            WelfareChange_HH_wealth(1,ww)=pctile.lb(ww);
            WelfareChange_HH_wealth(2,ww)=pctile.ub(ww);
            AggFVDist_HH = @(CEV) CompareVF_Indiv(V0, CEV, V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn_CEV, Params, Opt_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, SSvaluesFn, SSvalueParamNames, 1, pctile.lb(ww), pctile.ub(ww), ss, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter);
            try
            [WelfareChange_HH_wealth(3+ss,ww),CovergenceAlgorithms_HH_wealth(1+ss,ww)]=fminsearch(AggFVDist_HH,CEV_init,Options.fminoptions);
            disp(' ')
            fprintf('End of calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by wealth for income shock s = %.f \n', pctile.lb(ww), pctile.ub(ww), ss);
            if WelfareChange_HH_wealth(3+ss,ww) > 1
            CEV_init=WelfareChange_HH_wealth(3+ss,ww);
            elseif WelfareChange_HH_wealth(3+ss,ww) <= 1
            CEV_init=1.00;    
            end
            save Output/Reform/WelfareChange_HH_wealth.mat WelfareChange_HH_wealth CovergenceAlgorithms_HH_wealth
            xlswrite('Output\Tables\WelfareChange_HH_wealth.xlsx',WelfareChange_HH_wealth) % Export data for graphs in Excel
            catch
            WelfareChange_HH_wealth(3+ss,ww)=0;
            CovergenceAlgorithms_HH_wealth(1+ss,ww)=0;
            disp(' ')
            disp('CEV set to 0 since no agents in this part of the distribution');
            fprintf('End of calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by wealth for income shock s = %.f \n', pctile.lb(ww), pctile.ub(ww), ss);
            save Output/Reform/WelfareChange_HH_wealth.mat WelfareChange_HH_wealth CovergenceAlgorithms_HH_wealth
            xlswrite('Output\Tables\WelfareChange_HH_wealth.xlsx',WelfareChange_HH_wealth) % Export data for graphs in Excel
            end
        end
    end
end

% 7.2 Sort by income
n_pctile_groups = 10;
pctile.lb = (0:(100/n_pctile_groups):100-(100/n_pctile_groups))';
pctile.ub = pctile.lb+(100/n_pctile_groups);
pctile.ub(length(pctile.ub))=pctile.ub(length(pctile.ub))+0.00001; % To capture the very right tail of the distribution

if Options.Guess_HHtype == 1
    load Output/Reform/WelfareChange_HH_income.mat WelfareChange_HH_income CovergenceAlgorithms_HH_income
elseif Options.Guess_HHtype == 0
    WelfareChange_HH_income = ones(11,n_pctile_groups);
    CovergenceAlgorithms_HH_income = ones(9,n_pctile_groups);
end

SSvalueParamNames=struct();
SSvalueParamNames(1).Names={'e1','e2','e3','e4','w','r','J','omega'}; % Y
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J);
SSvaluesFn={SSvaluesFn_income};

for ww = 1:length(pctile.lb)
    for ss = 0:8 % 0 denotes all states
        if CovergenceAlgorithms_HH_income(1+ss,ww) > 10^(-3)
            IterCEVCounter=0;
            save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
            disp(' ')
            fprintf('Calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by income for income shock s = %.f... \n', pctile.lb(ww), pctile.ub(ww), ss);
            WelfareChange_HH_income(1,ww)=pctile.lb(ww);
            WelfareChange_HH_income(2,ww)=pctile.ub(ww);
            AggFVDist_HH = @(CEV) CompareVF_Indiv(V0, CEV, V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn_CEV, Params, Opt_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, SSvaluesFn, SSvalueParamNames, 1, pctile.lb(ww), pctile.ub(ww), ss, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter);
            try
            [WelfareChange_HH_income(3+ss,ww),CovergenceAlgorithms_HH_income(1+ss,ww)]=fminsearch(AggFVDist_HH,CEV_init,Options.fminoptions);
            disp(' ')
            fprintf('End of calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by income for income shock s = %.f \n', pctile.lb(ww), pctile.ub(ww), ss);
            if WelfareChange_HH_income(3+ss,ww) > 1
            CEV_init=WelfareChange_HH_income(3+ss,ww);
            elseif WelfareChange_HH_income(3+ss,ww) <= 1
            CEV_init=1.00;    
            end
            save Output/Reform/WelfareChange_HH_income.mat WelfareChange_HH_income CovergenceAlgorithms_HH_income
            xlswrite('Output\Tables\WelfareChange_HH_income.xlsx',WelfareChange_HH_income) % Export data for graphs in Excel
            catch
            WelfareChange_HH_income(3+ss,ww)=0;
            CovergenceAlgorithms_HH_income(1+ss,ww)=0;
            disp(' ')
            disp('CEV set to 0 since no agents in this part of the distribution');
            fprintf('End of calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by income for income shock s = %.f \n', pctile.lb(ww), pctile.ub(ww), ss);
            save Output/Reform/WelfareChange_HH_income.mat WelfareChange_HH_income CovergenceAlgorithms_HH_income
            xlswrite('Output\Tables\WelfareChange_HH_income.xlsx',WelfareChange_HH_income) % Export data for graphs in Excel
            end
        end
    end
end

% 7.3 Sort by VF
n_pctile_groups = 10;
pctile.lb = (0:(100/n_pctile_groups):100-(100/n_pctile_groups))';
pctile.ub = pctile.lb+(100/n_pctile_groups);
pctile.ub(length(pctile.ub))=pctile.ub(length(pctile.ub))+0.00001; % To capture the very right tail of the distribution

if Options.Guess_HHtype == 1
    load Output/Reform/WelfareChange_HH_VF.mat WelfareChange_HH_VF CovergenceAlgorithms_HH_VF
elseif Options.Guess_HHtype == 0
    WelfareChange_HH_VF = ones(11,n_pctile_groups);
    CovergenceAlgorithms_HH_VF = ones(9,n_pctile_groups);
end

% SSvalueParamNames and SSvaluesFn are taken from other definitions, here it does not matter

for ww = 1:length(pctile.lb)
    for ss = 0:8 % 0 denotes all states
        if CovergenceAlgorithms_HH_VF(1+ss,ww) > 10^(-3)
            IterCEVCounter=0;
            save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
            disp(' ')
            fprintf('Calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by Value Function for income shock s = %.f... \n', pctile.lb(ww), pctile.ub(ww), ss);
            WelfareChange_HH_VF(1,ww)=pctile.lb(ww);
            WelfareChange_HH_VF(2,ww)=pctile.ub(ww);
            AggFVDist_HH = @(CEV) CompareVF_Indiv(V0, CEV, V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn_CEV, Params, Opt_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, SSvaluesFn, SSvalueParamNames, 0, pctile.lb(ww), pctile.ub(ww), ss, Options.Tolerance, Options.Howards, Options.Howards2, Options.NCores, Options.Simperiods, Options.Burnin, Options.MaxIter);
            try
            [WelfareChange_HH_VF(3+ss,ww),CovergenceAlgorithms_HH_VF(1+ss,ww)]=fminsearch(AggFVDist_HH,CEV_init,Options.fminoptions);
            disp(' ')
            fprintf('End of calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by Value Function for income shock s = %.f \n', pctile.lb(ww), pctile.ub(ww), ss);
            if WelfareChange_HH_VF(3+ss,ww) > 1
            CEV_init=WelfareChange_HH_VF(3+ss,ww);
            elseif WelfareChange_HH_VF(3+ss,ww) <= 1
            CEV_init=1.00;    
            end
            save Output/Reform/WelfareChange_HH_VF.mat WelfareChange_HH_VF CovergenceAlgorithms_HH_VF
            xlswrite('Output\Tables\WelfareChange_HH_VF.xlsx',WelfareChange_HH_VF) % Export data for graphs in Excel
            catch
            WelfareChange_HH_VF(3+ss,ww)=0;
            CovergenceAlgorithms_HH_VF(1+ss,ww)=0;
            disp(' ')
            disp('CEV set to 0 since no agents in this part of the distribution');
            fprintf('End of calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f sorted by Value Function for income shock s = %.f \n', pctile.lb(ww), pctile.ub(ww), ss);
            save Output/Reform/WelfareChange_HH_VF.mat WelfareChange_HH_VF CovergenceAlgorithms_HH_VF
            xlswrite('Output\Tables\WelfareChange_HH_VF.xlsx',WelfareChange_HH_VF) % Export data for graphs in Excel
            end
        end
    end
end

elseif Options.SkipLoop_HHtype == 1
    disp(' ')
    disp('Loop to compute welfare changes (CEV) by type of household is not activated. Previous results are loaded.')
    load Output/Reform/WelfareChange_HH_VF.mat WelfareChange_HH_VF CovergenceAlgorithms_HH_VF
    load Output/Reform/WelfareChange_HH_wealth.mat WelfareChange_HH_wealth CovergenceAlgorithms_HH_wealth
    load Output/Reform/WelfareChange_HH_income.mat WelfareChange_HH_income CovergenceAlgorithms_HH_income
end

%% 8. Depict changes in distributions between baseline economy and optimal economy

% Share of the stationary distribution concentrated in each percentile group
% Note: we want to classify the households by wealth, by income, and by value function
n_pctile_groups = 10;
pctile.lb = (0:(100/n_pctile_groups):100-(100/n_pctile_groups))';
pctile.ub = pctile.lb+(100/n_pctile_groups);
pctile.ub(length(pctile.ub))=pctile.ub(length(pctile.ub))+0.00001; % To capture the very right tail of the distribution

% 8.1 Sort by wealth
StationaryDistrChange_HH_wealth = zeros(5,n_pctile_groups);

SSvalueParamNames=struct();
SSvalueParamNames(1).Names={};
SSvaluesFn_K = @(lab,anext,a,z) a; %K
SSvaluesFn={SSvaluesFn_K};

for ww = 1:length(pctile.lb)
    StationaryDistrChange_HH_wealth(1,ww)=pctile.lb(ww);
    StationaryDistrChange_HH_wealth(2,ww)=pctile.ub(ww);
    StationaryDistrChange_HH_wealth(3,ww)=GetStationaryDist_HHtype(V_BE, StationaryDistr_BE, Policy_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames, 1, (pctile.lb(ww)/100), (pctile.ub(ww)/100), 0);
    StationaryDistrChange_HH_wealth(4,ww)=GetStationaryDist_HHtype(V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params, SSvaluesFn, SSvalueParamNames, 1, (pctile.lb(ww)/100), (pctile.ub(ww)/100), 0);
    StationaryDistrChange_HH_wealth(5,ww)=(StationaryDistrChange_HH_wealth(4,ww)-StationaryDistrChange_HH_wealth(3,ww))/StationaryDistrChange_HH_wealth(3,ww)*100;
end

% 8.2 Sort by income
StationaryDistrChange_HH_income = zeros(5,n_pctile_groups);

SSvalueParamNames=struct();
SSvalueParamNames(1).Names={'e1','e2','e3','e4','w','r','J','omega'}; % Y
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J);
SSvaluesFn={SSvaluesFn_income};

for ww = 1:length(pctile.lb)
    StationaryDistrChange_HH_income(1,ww)=pctile.lb(ww);
    StationaryDistrChange_HH_income(2,ww)=pctile.ub(ww);
    StationaryDistrChange_HH_income(3,ww)=GetStationaryDist_HHtype(V_BE, StationaryDistr_BE, Policy_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames, 1, (pctile.lb(ww)/100), (pctile.ub(ww)/100), 0);
    StationaryDistrChange_HH_income(4,ww)=GetStationaryDist_HHtype(V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params, SSvaluesFn, SSvalueParamNames, 1, (pctile.lb(ww)/100), (pctile.ub(ww)/100), 0);
    StationaryDistrChange_HH_income(5,ww)=(StationaryDistrChange_HH_income(4,ww)-StationaryDistrChange_HH_income(3,ww))/StationaryDistrChange_HH_income(3,ww)*100;
end

% 8.3 Sort by VF
StationaryDistrChange_HH_VF = zeros(5,n_pctile_groups);

% SSvalueParamNames and SSvaluesFn are taken from other definitions, here it does not matter

for ww = 1:length(pctile.lb)
    StationaryDistrChange_HH_VF(1,ww)=pctile.lb(ww);
    StationaryDistrChange_HH_VF(2,ww)=pctile.ub(ww);
    StationaryDistrChange_HH_VF(3,ww)=GetStationaryDist_HHtype(V_BE, StationaryDistr_BE, Policy_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames, 0, (pctile.lb(ww)/100), (pctile.ub(ww)/100), 0);
    StationaryDistrChange_HH_VF(4,ww)=GetStationaryDist_HHtype(V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params, SSvaluesFn, SSvalueParamNames, 0, (pctile.lb(ww)/100), (pctile.ub(ww)/100), 0);
    StationaryDistrChange_HH_VF(5,ww)=(StationaryDistrChange_HH_VF(4,ww)-StationaryDistrChange_HH_VF(3,ww))/StationaryDistrChange_HH_VF(3,ww)*100;
end

% 8.4 Save results on disk
save Output/Reform/StationaryDistrChange_HH.mat StationaryDistrChange_HH_wealth StationaryDistrChange_HH_income StationaryDistrChange_HH_VF

%% 9. Changes in hours worked (between baseline and optimal economy), ranked by earnings

% Average hours worked in each percentile group
% Note: we want to classify the households by earnings
n_pctile_groups = 10;
pctile.lb = (0:(100/n_pctile_groups):100-(100/n_pctile_groups))';
pctile.ub = pctile.lb+(100/n_pctile_groups);
pctile.ub(length(pctile.ub))=pctile.ub(length(pctile.ub))+0.00001; % To capture the very right tail of the distribution

HoursChange = zeros(5,n_pctile_groups);

SSvalueParamNames(1).Names={};
SSvaluesFn_H = @(lab,anext,a,z) lab; %H
SSvalueParamNames(2).Names={'e1','e2','e3','e4','w','r','J'};
SSvaluesFn_earnings = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r; % Earnings
SSvaluesFn={SSvaluesFn_H,SSvaluesFn_earnings};

for ww = 1:length(pctile.lb)
    HoursChange(1,ww)=pctile.lb(ww);
    HoursChange(2,ww)=pctile.ub(ww);
    Outcome =  SSvalues_Indiv_Vars_Ev(V_BE, StationaryDistr_BE, Policy_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames);
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest (earnings in this case)
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=(pctile.lb(ww)/100)).*(Outcome(:,size(Outcome,2))<(pctile.ub(ww)/100))]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    HoursChange(3,ww)=sum(sum((Outcome(:,5)).*Outcome(:,size(Outcome,2)))); % Compute aggregate value function for individuals in the percentile of interestShare = sum(sum(Outcome(:,size(Outcome,2))));
    Outcome =  SSvalues_Indiv_Vars_Ev(V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params, SSvaluesFn, SSvalueParamNames);
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest (earnings in this case)
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=(pctile.lb(ww)/100)).*(Outcome(:,size(Outcome,2))<(pctile.ub(ww)/100))]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    HoursChange(4,ww)=sum(sum((Outcome(:,5)).*Outcome(:,size(Outcome,2)))); % Compute aggregate value function for individuals in the percentile of interestShare = sum(sum(Outcome(:,size(Outcome,2))));
    HoursChange(5,ww)=(HoursChange(4,ww)-HoursChange(3,ww))/HoursChange(3,ww)*100;
end

% Save results on disk
save Output/Reform/HoursChange.mat HoursChange

%% 10. Changes in savings rates (between baseline and optimal economy), ranked by wealth

% Average savings rate in each percentile group
% Note: we want to classify the households by wealth
n_pctile_groups = 10;
pctile.lb = (0:(100/n_pctile_groups):100-(100/n_pctile_groups))';
pctile.ub = pctile.lb+(100/n_pctile_groups);
pctile.ub(length(pctile.ub))=pctile.ub(length(pctile.ub))+0.00001; % To capture the very right tail of the distribution

SavingsChange = zeros(5,n_pctile_groups);

SSvalueParamNames(1).Names={'delta'};
SSvaluesFn_I = @(lab,anext,a,z,delta) anext-a*(1-delta); %I
SSvalueParamNames(2).Names={'e1','e2','e3','e4','w','r','J','omega'}; 
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J); % income
SSvalueParamNames(3).Names={};
SSvaluesFn_K = @(lab,anext,a,z) a; %K
SSvaluesFn={SSvaluesFn_I,SSvaluesFn_income,SSvaluesFn_K};

for ww = 1:length(pctile.lb)
    SavingsChange(1,ww)=pctile.lb(ww);
    SavingsChange(2,ww)=pctile.ub(ww);
    Outcome =  SSvalues_Indiv_Vars_Ev(V_BE, StationaryDistr_BE, Policy_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames);
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest (wealth in this case)
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=(pctile.lb(ww)/100)).*(Outcome(:,size(Outcome,2))<(pctile.ub(ww)/100))]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    Outcome=[Outcome (Outcome(:,5))./(Outcome(:,6))];
    Outcome(isnan(Outcome))=0;
    Outcome=[Outcome Outcome(:,size(Outcome,2)).*Outcome(:,size(Outcome,2)-1)];
    SavingsChange(3,ww)=sum(Outcome(:,size(Outcome,2)));
    Outcome =  SSvalues_Indiv_Vars_Ev(V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params, SSvaluesFn, SSvalueParamNames);
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest (wealth in this case)
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=(pctile.lb(ww)/100)).*(Outcome(:,size(Outcome,2))<(pctile.ub(ww)/100))]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    Outcome=[Outcome (Outcome(:,5))./(Outcome(:,6))];
    Outcome(isnan(Outcome))=0;
    Outcome=[Outcome Outcome(:,size(Outcome,2)).*Outcome(:,size(Outcome,2)-1)];
    SavingsChange(4,ww)=sum(Outcome(:,size(Outcome,2)));
    SavingsChange(5,ww)=(SavingsChange(4,ww)-SavingsChange(3,ww))/SavingsChange(3,ww)*100;
end

% Save results on disk
save Output/Reform/SavingsChange.mat SavingsChange

%% 11. Changes in average effective tax rate (between baseline and optimal economy), ranked by income

% Average effective tax rate in each percentile group
% Note: we want to classify the households by wealth
pctile.lb=[1  1  21 31 21 41 41 51 61 61 71 81 81     90 95 99    ];
pctile.ub=[10 20 30 40 40 60 50 60 80 70 80 90 100.01 95 99 100.01]; % 100.01 o capture the very right tail of the distribution

EffTaxRateChange = zeros(5,length(pctile.lb));

SSvalueParamNames(1).Names={'J','r','w','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
SSvaluesFn_tax = @(lab,anext,a,z,J,r,w,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) (TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa));
SSvalueParamNames(2).Names={'e1','e2','e3','e4','w','r','J','omega'}; % income
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J);
SSvaluesFn={SSvaluesFn_tax,SSvaluesFn_income};

for ww = 1:length(pctile.lb)
    EffTaxRateChange(1,ww)=pctile.lb(ww);
    EffTaxRateChange(2,ww)=pctile.ub(ww);
    Outcome =  SSvalues_Indiv_Vars_Ev(V_BE, StationaryDistr_BE, Policy_BE, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames);
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest (wealth in this case)
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=(pctile.lb(ww)/100)).*(Outcome(:,size(Outcome,2))<(pctile.ub(ww)/100))]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    Outcome=[Outcome (Outcome(:,5))./(Outcome(:,6))];
    Outcome(isnan(Outcome))=0;
    Outcome=[Outcome Outcome(:,size(Outcome,2)).*Outcome(:,size(Outcome,2)-1)];
    EffTaxRateChange(3,ww)=sum(Outcome(:,size(Outcome,2)));
    Outcome =  SSvalues_Indiv_Vars_Ev(V_opt, StationaryDistr_opt, Policy_opt, n_d, n_a, n_z, d_grid, a_grid, z_grid, Opt_Params, SSvaluesFn, SSvalueParamNames);
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest (wealth in this case)
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=(pctile.lb(ww)/100)).*(Outcome(:,size(Outcome,2))<(pctile.ub(ww)/100))]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    Outcome=[Outcome (Outcome(:,5))./(Outcome(:,6))];
    Outcome(isnan(Outcome))=0;
    Outcome=[Outcome Outcome(:,size(Outcome,2)).*Outcome(:,size(Outcome,2)-1)];
    EffTaxRateChange(4,ww)=sum(Outcome(:,size(Outcome,2)));
    EffTaxRateChange(5,ww)=(EffTaxRateChange(4,ww)-EffTaxRateChange(3,ww))/EffTaxRateChange(3,ww);
end

% Save results on disk
save Output/Reform/EffTaxRateChange.mat EffTaxRateChange
xlswrite('Output\Tables\EffTaxRateChange.xlsx',EffTaxRateChange) % Export data for graphs in Excel
