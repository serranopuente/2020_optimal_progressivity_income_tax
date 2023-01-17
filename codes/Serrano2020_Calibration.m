% Optimal Progressivity of Personal Income Tax: A General Equilibrium Evaluation for Spain
% Calibration Algorithm
% Darío Serrano-Puente
% Banco de España | DG Economics, Statistics and Research | Structural Analysis Division
% darioserrapuente@gmail.com | dario.serrano@bde.es
% https://sites.google.com/view/darioserranopuente/
% September 10th, 2020
%**********************************************************************
% Name: Main file with initial parametrization and calibration
% File name: Serrano2020_Calibration.m
% Connected files: in subfolder named "subroutines"

% Description: solution algorithm for the modfified version of the 
% stochastic neoclassical growth model with uninsured idiosyincratic risk
% and no aggregate risk, inspired on the Aiyagari (1994,QJE) model. Based
% on the code originally written by Robert Kirkby to replicate the paper
% by Castaneda, Diaz-Gimenez, and Rios-Rull (2003).
%**********************************************************************
clc;            % Clear screen
clear;          % Clear memory
rng default;    % Fix the seed for the random number generator (so as to keep simulations the same)
close all;      % Close open windows

addpath(genpath(pwd));
%**********************************************************************

%% 1. Set some options for the code

disp('Running "Optimal Progressivity of Personal Income Tax: A General Equilibrium Evaluation for Spain"')
disp('Calibration Algorithm')
disp('Serrano-Puente, Darío (2020)')
disp('Current code version: September 10th, 2020')
disp(' ')
% Parallelize on all possible cores of the computer
delete(gcp('nocreate')); % Close existing interactive parpool session if so
parpool; % Use 2 cores (number of available cores of my computer)
% parpool(1) % Just use 1 core if this line is commented out
PoolDetails=gcp; % Get the number of cores working
Options.SkipGE=1; % Just a placeholder to work on codes without rerunning the baseline General Equilibrium step - 0:='run the GE finder algorithm'; 1:='rertieve GE results from previous finding'
Options.MSM=0; % Placeholder to work on codes by 'running the Method of Simulated Moments calibration procedure' (1) or 'retrieving last calibration results / manual calibration' (0)
Options.Test=0; % To test whether some functions work properly
Options.NCores=PoolDetails.NumWorkers; % Number of workers/cores to use in the stationary distribution simulation
Options.Tolerance=10^(-5); % Iterate until reach this convergence level
Options.Howards=60;
Options.Howards2=500;
Options.Burnin=1000;
Options.Simperiods=10^5; % Simulation periods to get the stationary distribution
Options.NSim=8*10^5; % Used when calculating the statistics on intergenerational earnings correlation and on life-cycle earnings profile. If less than 10^5~10^6, the moments bounce around too much (not enough simulations to really pin it down).
Options.MaxIter=40000;
n_l=13;  % Number of grid points for labor supply
% Two General Equilibrium variables: interest rate (r), linear tax rate of remaining taxes (kappa)
n_r=551; % Number of grid points to search for the equilibirum interest rate - Irrelevant for this type of GE-finding algorithm
n_kappa=41; % Number of grid points to search for tax rate in equilibirum  - Irrelevant for this type of GE-finding algorithm
            % Note: this parameter 'kappa' is properly calibrated when finding GE, not just when calibration (MSM) is done.
disp(' ')
disp('Code options:')
fprintf('Skip General Equilibrium calculation: %.f \n', Options.SkipGE)
fprintf('Running calibration algorithm: %.f \n', Options.MSM)
fprintf('Number of cores connected to: %.f \n', Options.NCores)
fprintf('Tolerance level: %.5f \n', Options.Tolerance)
fprintf('Distance of Howards Value Function Iteration Improvement: %.f \n', Options.Howards2)
fprintf('Maximum number of iterations of Howards Value Function Iteration Improvement: %.f \n', Options.Howards)
fprintf('Burnin: %.f \n', Options.Burnin)
fprintf('Simulation periods to get the stationary distribution: %.f \n', Options.Simperiods)
fprintf('Number of simulations for intergenerational correlation and life-cycle profile: %.f \n', Options.NSim)
fprintf('Maximum number of iterations: %.f \n', Options.MaxIter)

%% 2. Calibration targets (statistics from micro-data)
% All statistics depict a figure of the Spanish economy in 2015

% 2.0 Targets to match (22 targets)
CalibStats_Data=zeros(22,1);

EstimationTargetNames={'Capital-output ratio','Gov. expenditure to output Ratio','Transfers to output ratio','Share of disposable time allocated to market','Ratio of earnings old to young','Correlation of income between fathers and sons',...
   'Income Gini', 'Wealth Gini','Income share up to p40','Income share p41-p60','Income share p61-p80','Income share p81-p100','Wealth share up to p40','Wealth share p41-p60','Wealth share p61-p80','Wealth share p81-p100',...
   'Income share p90-p95','Income share p95-p99','Income share p99-p100','Wealth share p90-p95','Wealth share p95-p99','Wealth share p99-p100'};

% 2.1 Macroeconomic aggregates
EstimationTargets.CapitalOutputRatio=4.2508; % K/Y
CalibStats_Data(1)=EstimationTargets.CapitalOutputRatio;
EstimationTargets.CapitalIncomeShare=0.4755;
% Params.alpha=0.4755; % Follows immediately from labor income share.
EstimationTargets.InvestmentToOutputRatio=0.2194; % I/Y
% Total tax revenues as percentage of GDP in Spain in 2015 was 0.3363.
EstimationTargets.GovExpenditureToOutputRatio=100*0.2227; % G/Y
CalibStats_Data(2)=EstimationTargets.GovExpenditureToOutputRatio;
EstimationTargets.TransfersToOutputRatio=100*0.1136; %Tr/Y
CalibStats_Data(3)=EstimationTargets.TransfersToOutputRatio;

% 2.2 Allocation of time and consumption
% Params.elle=3.2;
EstimationTargets.ShareOfDisposableTimeAllocatedToMarket=100*0.3083; % Avg. share of disposable time allocated to working in the market
CalibStats_Data(4)=EstimationTargets.ShareOfDisposableTimeAllocatedToMarket;

% 2.3 Age structure of the population
EstimationTargets.ExpectedDurationOfWorkingLife=35.0000;
EstimationTargets.ExpectedDurationOfRetirement=22.7661;
% These lead us directly to
% Params.p_eg=0.0286; % Note: 1/p_eg=35.0000
% Params.p_gg=0.9561; % Note: 1/(1-p_gg)=22.7661
% Follows from theoretical results on 'survival analysis': the expected duration of process with constant-hazard-rate lambda is 1/lambda. 
% Here p_eg and (1-p_gg) are the hazard rates. See, e.g., example on middle of pg 3 of http://data.princeton.edu/wws509/notes/c7.pdf.

% 2.4 Life-cycle profile of earnings
EstimationTargets.RatioOfEarningsOldtoYoung=1.5620; % Ratio of average earnings for households between ages of 41 & 65 to average earnings of households between ages of 18 & 40.
CalibStats_Data(5)=EstimationTargets.RatioOfEarningsOldtoYoung;

% 2.5 Intergenerational transmission of earnings ability
EstimationTargets.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=0.5044; % Cross-sectional correlation between avg. life-time earnings of one generation of HHs and the avg. life-time earnings of their immediate descendants
CalibStats_Data(6)=EstimationTargets.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons;

% 2.6 Income taxation
EstimationTargets.GovernmentBudgetBalance=0; % This is a General Equilibrium condition, so no need to repeat it here.

% 2.7 Normalization
% Params.e1=1; % One is actually better off normalizing Params.e2=1 than Params.e1=1.
% Normalize the diagonal elements of Gamma_ee (ie., Gamma_ee_11,
% Gamma_ee_22, Gamma_ee_33, Gamma_ee_44). Choose these as simply setting,
% e.g., Gamma_ee_11=1-Gamma_ee_12-Gamma_ee_13-Gamma_ee_14-p_eg is unlikely to
% lead us to a negative value of Gamma_ee_11. Thus in terms of maintaining
% the constraints on Gamma_ee (all elements between zero and one, rows sum
% to one) required for it to be a transition matrix, we are less likely to
% be throwing out parameter vectors when estimating because they failed to
% meet these constraints. This is just a numerical trick that works well in
% practice as we 'know' that most of the weight of the transition matrix is
% on the diagonal.
% Note that this normalization of the diagonal elements of Gamma_ee is
% actually hard-coded into the how we have written the codes that create
% the transition matrix.

% 2.8 Distributions of Income and wealth
EstimationTargets.IncomeGini=0.4800;
CalibStats_Data(7)=EstimationTargets.IncomeGini;
EstimationTargets.WealthGini=0.6820;
CalibStats_Data(8)=EstimationTargets.WealthGini;
EstimationTargets.IncomeQuintileSharesAsFraction=100*[0.1272, 0.1384, 0.2119, 0.5225]; % Percentiles shares: up to p40, p40-p60, p60-p80, p80-100
CalibStats_Data(9:12)=EstimationTargets.IncomeQuintileSharesAsFraction;
EstimationTargets.WealthQuintileSharesAsFraction=100*[0.0362, 0.0965, 0.1811, 0.6862]; % Percentiles shares: up to p40, p40-p60, p60-p80, p80-100
CalibStats_Data(13:16)=EstimationTargets.WealthQuintileSharesAsFraction;
EstimationTargets.IncomeTopSharesAsFraction=100*[0.1101,0.1340,0.1207]; % Shares of p90-p95, p95-p99, p99-p100
CalibStats_Data(17:19)=EstimationTargets.IncomeTopSharesAsFraction;
EstimationTargets.WealthTopSharesAsFraction=100*[0.1293,0.1979,0.2027]; % Shares of p90-p95, p95-p99, p99-p100
CalibStats_Data(20:22)=EstimationTargets.WealthTopSharesAsFraction;

% 2.9 Save the calibration targets
save ./Output/Calib/CalibStats_Data.mat CalibStats_Data
NumOfCalibStats=length(CalibStats_Data); % Number of calibration targets
disp(' ')
fprintf('Number of moments/statistics to be matched by the calibration procedure: %.f \n', NumOfCalibStats)

%% 3. Parametrization

% 3.0 Number of substates when retired or working-age agent
Params.J=4; % The code is just developed to support J = 4.

% 3.1 Preferences and utility
Params.beta=0.95510;   % Subjective time discount factor
Params.gamma=1.5000;   % Curvature of consumption - Inverse elasticity of intertemporal substitution as standard in literature - CRRA (must be != 1)
Params.varphi=2.6500;  % Curvature of leisure - Frisch elasticity of labor supply as standard in literature (multiplicative inverse of gamma)
Params.chi=0.500164;   % Relative share of consumption and leisure
Params.elle=3.2;       % Productive time

% 3.2 Age and employment process
Params.p_eg=0.0286;    % Common probability of retiring (Eurostat) - 35.1000 is expected duration of working life
Params.p_gg=0.9561;    % 1 - p_gg = Common probability of dying - 23.4527 is expected duration of retirement
Params.phi1=0.9999;   % Earnings life cycle controller
Params.phi2=0.9715;    % Intergenerational earnings persistence controller

% 3.3 Technology
Params.alpha=0.4755;   % Capital income share - (1-labor_share)
Params.delta=0.0516;   % Depreciation rate of physical capital - Follows immediately from delta=I/K in stationary general eqm, hence delta=(I/Y)/(K/Y)

% 3.4 Government Policy
Params.G=0.2227;       % Government expenditures. Note: G is not really a true 'parameter' to be calibrated, but just an auxiliary parameter that represtents G/Y in the GE calculation and is intended to acquire government budget balance which is a 'market clearance' condition - Then, after the calibration and the GE calculation it it will be replaced by G or G/Y depending on what magnitude the planner wants to keep unchanged in the reformed scenario
Params.omega=3.21500;  % Normalized transfers to retirees
Params.lambda=0.8924;  % Average taxes - HSV (2017) specification - Estimation à la Roberto Ramos (2019)
Params.tau=0.1146;     % Progressivity of PIT system - HSV (2017) specification - Estimation à la Roberto Ramos (2019)
Params.kappa=0.0310;   % Linear term of remaining taxes

% 3.5 Process on exogenous shocks
Params.e1=1.00; Params.e2=2.708056; Params.e3=7.800000; Params.e4=90.000000;
Params.e=[Params.e1,Params.e2,Params.e3,Params.e4,0,0,0,0];
% In Castañeda et al. (2003), they do the nomalization on e(2). So e1=1/e2; e3=e3/e2; e4=e4/e2; e2=1; % No longer need to do this.
% Note: the diagonal elements of Gamma_ee can be considered as residually determined to make each row of Gamma (Gamma_ee) add up to 1 (to 1-p_eg); 
% In principle it doesn't matter which element in each row is residual, but from perspective of calibration it is easiest to 
% use the diagonals to avoid contraint that they as transition probabilities they must be >=0 from binding.
Params.Gamma_ee_12=0.100605; Params.Gamma_ee_13=0.000100; Params.Gamma_ee_14=0.000500;     % Params.Gamma_ee_11=~0.95; 
Params.Gamma_ee_21=0.023541; Params.Gamma_ee_23=0.010000; Params.Gamma_ee_24=0.000100;     % Params.Gamma_ee_22=~0.97;
Params.Gamma_ee_31=0.000100; Params.Gamma_ee_32=0.035000; Params.Gamma_ee_34=0.000419;     % Params.Gamma_ee_33=~0.98;
Params.Gamma_ee_41=0.000100; Params.Gamma_ee_42=0.016819; Params.Gamma_ee_43=0.000100;     % Params.Gamma_ee_44=~0.75;

% 3.6 Determine parameters directly identified (they will not be estimated/calibrated)
DirectParamNames={'J','gamma','elle','p_eg','p_gg','alpha','delta','G','lambda','tau','e1'};
    for ii=1:length(DirectParamNames)
        DirectParams.(DirectParamNames{ii})=Params.(DirectParamNames{ii});
    end

% 3.7 Determine parameters to estimate (22 parameters) - the rest will be direct parameters (they will not be estimated as they will take given values)
ParamNamesToEstimate={'beta','varphi','chi','phi1','phi2','omega','kappa','e2','e3','e4',...
   'Gamma_ee_12','Gamma_ee_13','Gamma_ee_14','Gamma_ee_21','Gamma_ee_23','Gamma_ee_24','Gamma_ee_31','Gamma_ee_32','Gamma_ee_34','Gamma_ee_41','Gamma_ee_42','Gamma_ee_43'};
% Note: 'r' and 'kappa' are determined by the General Equilibrium conditions, rather than the calibration.
NumOfParamsToEstimate=length(ParamNamesToEstimate);
disp(' ')
fprintf('Number of model parameters to be calibrated: %.f \n', NumOfParamsToEstimate)

% 3.8 Bounds on parameter values (search the estimated parameters must be inside these intervals)
ParamBounds.lb=[0.95505,   2.6400,  0.5001,  0.9998, 0.9700,  3.2050,     0.0300,  2.6500,  7.7000,   90.0000,          0.0950,         0.0001,         0.0001,         0.0200,         0.0050,         0.0001,         0.0001,         0.0250,         0.0001,         0.0001,         0.0140,         0.0001]';
ParamBounds.ub=[0.95515,   2.6600,  0.5002,  0.9999, 0.9730,  3.2200,     0.0325,  2.7500,  8.0000,  100.0000,          0.1100,         0.0002,         0.0008,         0.0280,         0.0120,         0.0002,         0.0002,         0.0360,         0.0005,         0.0002,         0.0200,         0.0002]';
             % ['beta', 'varphi',   'chi',  'phi1', 'phi2', 'omega',    'kappa',    'e2',    'e3',      'e4',   'Gamma_ee_12',  'Gamma_ee_13',  'Gamma_ee_14',  'Gamma_ee_21',  'Gamma_ee_23',  'Gamma_ee_24',  'Gamma_ee_31',  'Gamma_ee_32',  'Gamma_ee_34',  'Gamma_ee_41',  'Gamma_ee_42',  'Gamma_ee_43']
% Note: lb & ub are both column vectors.

% 3.9 Type of distance between estimated moment and targeted moment to minimize (1:='squared difference'; 2:='absolute difference')
CalibDistType=ones(NumOfCalibStats,1)';
estimationoptions.TargetDistanceFns.IncomeQuintileSharesAsFraction=1;
CalibDistType(9:12)=estimationoptions.TargetDistanceFns.IncomeQuintileSharesAsFraction;
estimationoptions.TargetDistanceFns.WealthQuintileSharesAsFraction=1;
CalibDistType(13:16)=estimationoptions.TargetDistanceFns.WealthQuintileSharesAsFraction;

% 3.10 Calibration weights - Relative importance of the parameters in the calibration process
CalibWeights=ones(NumOfCalibStats,1);
TargetWeights.CapitalIncomeRatio=1000;
CalibWeights(1)=TargetWeights.CapitalIncomeRatio;
TargetWeights.TransfersToOutputRatio=20;
CalibWeights(3)=TargetWeights.TransfersToOutputRatio;
TargetWeights.GovernmentExpToOutputRatio=20;
CalibWeights(2)=TargetWeights.GovernmentExpToOutputRatio;
% Targets include an excess of inequality statistics, so the weights given to these do not need to be very high.
% However, they must be increased a bit as they are important part of the purpose of model, and were being ignored during the calibration otherwise (in earlier runs).
TargetWeights.IncomeQuintileSharesAsFraction=15;
CalibWeights(9:12)=TargetWeights.IncomeQuintileSharesAsFraction*ones(length(TargetWeights.IncomeQuintileSharesAsFraction),1);
TargetWeights.IncomeTopSharesAsFraction=30;
CalibWeights(17:19)=TargetWeights.IncomeTopSharesAsFraction*ones(length(TargetWeights.IncomeTopSharesAsFraction),1);
TargetWeights.WealthQuintileSharesAsFraction=15;
CalibWeights(13:16)=TargetWeights.WealthQuintileSharesAsFraction*ones(length(TargetWeights.WealthQuintileSharesAsFraction),1);
TargetWeights.WealthTopSharesAsFraction=30;
CalibWeights(20:22)=TargetWeights.WealthTopSharesAsFraction*ones(length(TargetWeights.WealthTopSharesAsFraction),1);
% The data and link to model are not strongest for the following two, so I give them lower weights.
TargetWeights.RatioOfEarningsOldtoYoung=40;
CalibWeights(5)=TargetWeights.RatioOfEarningsOldtoYoung;
TargetWeights.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=40;
CalibWeights(6)=TargetWeights.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons;
% An early estimation attempt ended up going off-track and making almost nobody work. Following makes the fraction of time worked estimation target important.
TargetWeights.ShareOfDisposableTimeAllocatedToMarket=20;
CalibWeights(4)=TargetWeights.ShareOfDisposableTimeAllocatedToMarket;         

%% 4. Grids

% 4.1 Markov-chain transition probability matrix, invariant probability measures and endowment grid/vector
[e_grid,Gamma,gammastar_initial,gammastarfull_initial]=Create_Exog_Shock(Params);

% 4.2 Labor supply grid
l_grid=linspace(0,Params.elle,n_l)';

% 4.3 Capital grid
% Proposed grid
a_grid=[0:0.025:1,1.0625:0.0625:2,2.16:0.16:10,10.6:0.6:100,105:4:1500]';
% A more accurate grid
% a_grid=[0:0.02:1, 1.05:0.05:2, 2.1:0.1:50, 50.5:0.5:100, 104:4:2000]';
% Values they used in Castañeda et al. (2003)
% a_grid=[0:0.02:1,1.05:0.05:2,2.1:0.1:10,10.5:0.5:100,104:4:1500]';

% 4.4 Grid of the two endogenous variables together
d_grid=[l_grid; a_grid]; % Notational conventions used by VFI subroutines

% 4.5 Grid of shock number (determines age & retirement)
z_grid=linspace(1,2*Params.J,2*Params.J)';

% 4.6 Display sizes to get an idea of the dimensions of the Dynamic Programming problem
n_a=length(a_grid); % Number of grid points for assets/capital
n_d=[n_l,n_a];      % Length of endogenous grids
n_z=length(z_grid); % Length of state grid
n_p=[n_r,n_kappa];  % Length of solvers' grid
disp(' ')
disp('Grid sizes:')
fprintf('Number of grid points for labor choice: %.f \n', n_d(1))
fprintf('Number of grid points for assets/capital: %.f \n', n_a)
fprintf('Length of state (age) grid: %.f \n', n_z)
N_d=prod(n_d); N_a=prod(n_a); N_s=prod(n_z);

% 4.7 Save the grids
save ./Output/Grids.mat d_grid a_grid z_grid n_d n_a n_z

%% 5. Utility function and next period's assets from current period's decisions

% 5.1 Define utility (or return) function equation by calling the return's calculation subroutine
UtilityFn=@(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa) ReturnFn(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa);
UtilityFnParamNames={'r','gamma','varphi','chi','elle','alpha','delta','e1','e2','e3','e4','omega','lambda','tau','kappa'};
DiscountFactorParamNames={'beta'};

% 5.2 Determine next period's assets from current period's decisions
PhiaprimeParamNames={};
Phi_aprimeMatrix=PhiaprimeMatrix(n_d,n_z);

%% 6. Pre-calibration

disp(' ')
disp('Initializing pre-calibration...')

% 6.1 Set interest rate to be what must be the eqm value ('Handy Trick' to avoid need to loop over the calculation of GE during calibration process)
disp(' ')
disp('Setting interest rate in equilibirium conditions directily from the data.')
p_eqm(1)=Params.alpha*1/CalibStats_Data(1)-Params.delta;
Params.r=p_eqm(1);
Params.w=(1-Params.alpha)*(((Params.r+Params.delta)/(Params.alpha))^(Params.alpha/(Params.alpha-1)));

% 6.2 Get the VF and the stationary distribution from the initial guess
disp(' ')
disp('Retrieving value function and stationary dsitribution from initial guess...')
V0=ones(n_a,n_z); % Initial guess for the value function
tic;
[V, Policy]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2);
rtimes(1)=toc;
tic;
StationaryDistr=StationaryDist(Policy,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,Options.NCores,Options.Simperiods,Options.Burnin,Options.Tolerance,Options.MaxIter);
rtimes(2)=toc;
    
% 6.3 Get the model statistics
tic;
ModelStats_initvals=ModelStatistics(Options.NSim,StationaryDistr, Policy, Phi_aprimeMatrix, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Params);
rtimes(3)=toc;
ExtraStats_initvals=ModelStatistics_ExtraStats(StationaryDistr, Policy, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params);
    
% 6.4 Save pre-calibration results
for ii=1:NumOfParamsToEstimate
    CalibParams_initvals.(ParamNamesToEstimate{ii})=Params.(ParamNamesToEstimate{ii});
end
CalibParams_initvalsVec=CreateVectorFromParams(Params, ParamNamesToEstimate)';
save Output/Calib/Pre-Calib_Results.mat p_eqm V Policy a_grid StationaryDistr DirectParams CalibParams_initvals CalibParams_initvalsVec ModelStats_initvals ExtraStats_initvals rtimes
disp(' ')
fprintf('Run times: VFI: %.4f, Stationary Distr.: %.4f, Baseline Model Statistics: %.4f \n',rtimes(1), rtimes(2), rtimes(3))

% 6.5 Check interest rates
disp(' ')
fprintf('So we have p_eqm (calculated from K/Y from data) as %.4f \n', p_eqm(1))  % alpha*1/CalibStats_Data(1)-delta
SSvalueParamNames(1).Names={};
SSvaluesFn_1 = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
SSvaluesFn_2 = @(lab,anext,a,z,e1,e2,e3,e4) lab*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvaluesFn={SSvaluesFn_1, SSvaluesFn_2};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDistr, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);
p_market=Params.alpha*(SSvalues_AggVars(1)^(Params.alpha-1))*(SSvalues_AggVars(2)^(1-Params.alpha))-Params.delta;
fprintf('While p_market (calculated from marginal product of capital) as %.4f \n', p_market)
Y=(SSvalues_AggVars(1)^Params.alpha)*(SSvalues_AggVars(2)^(1-Params.alpha));
K=SSvalues_AggVars(1);
L=SSvalues_AggVars(2);
fprintf('This is because: \n')
fprintf('Y = %.4f, K = %.4f, L = %.4f \n', Y,K,L)
fprintf('alpha = %.4f, delta = %.4f \n', Params.alpha, Params.delta)
p_eqm_check=Params.alpha*1/(K/Y)-Params.delta;
fprintf('Double check of p_market (based on K/Y): %.4f \n', p_eqm_check)
p_eqm_check=Params.alpha*1/(ModelStats_initvals(1))-Params.delta;
fprintf('Double check of p_market (based on Model Statistics): %.4f \n', p_eqm_check)
fprintf('Note: Both the checks of p_market should give the same value as p_market \n')

disp(' ')
fprintf('End of pre-calibration. \n')

%% 7. Calibration procedure: Minimizing the distance to the moments (Method of Simulated Moments, MSM)

disp(' ')
fprintf('Initializing calibration... \n')
CalibCounter=0; % Create and save a counter, then inside calibration loop one can load, increment, and resave, thus keeping count.
save ./Output/Calib/CalibCounter.mat CalibCounter

% 7.1 Set up function that takes parameters in and gives model statistics out
CalibStats_Model_Fn =@(CalibParamsVec) ModelTargetsFn(CalibParamsVec,Params,n_d,n_a,n_z,a_grid,UtilityFn,UtilityFnParamNames, DiscountFactorParamNames,PhiaprimeParamNames,ParamNamesToEstimate,NumOfCalibStats,EstimationTargetNames,CalibStats_Data,Options.NCores,Options.NSim,Options.Howards,Options.Howards2,Options.Simperiods,Options.Burnin,Options.Tolerance,Options.MaxIter);
if Options.Test==1
    disp(' ')
    disp('Test Model Targets Function')
    CalibStatsTest=CalibStats_Model_Fn(CalibParams_initvalsVec); %This is just to test it
    disp(' ')
    disp('End of Test Model Targets Function')
    CalibCounter=0; % Reset counter to start from 0 the calibration loop
    save ./Output/Calib/CalibCounter.mat CalibCounter
end

% 7.2 Define the function to be minimized: the sum of the weighted distances between model and data moments
CalibFn=@(CalibParamsVec) sum(MethodSimMomentsFn(CalibParamsVec, CalibStats_Model_Fn, CalibStats_Data, CalibDistType, CalibWeights, NumOfCalibStats));
if Options.Test==1
    disp(' ')
    disp('Test Calibration Function')
    DistTestFull=CalibFn(CalibParams_initvalsVec);
    disp(' ')
    disp('End of Test Calibration Function')
    CalibCounter=0; % Reset counter to start from 0 the calibration loop
    save ./Output/Calib/CalibCounter.mat CalibCounter
end

% 7.3 Prepare the settings for the CMA-ES algorithm, run it (if option activated) and save the final calibrated parameters
if Options.MSM==1
    
    % 7.3.1 Calibration loop (in case it is activated)
    disp(' ')
    disp('Starting calibration loop: CMA-ES algorithm of Andreasen (2010) [Covariance Matrix Adaptation Evolutionary Stategy]')
    tic;
    CMAES_Insigma=0.5*abs(ParamBounds.ub-ParamBounds.lb); % The std deviation in the initial search distributions
    CMAES_sigma=1;                                        % The step size
    opts.SigmaMax=5;                                      % The maximal value for sigma
    opts.LBounds=ParamBounds.lb;                          % Lower bound for CalibParams
    opts.UBounds=ParamBounds.ub;                          % Upper bound for CalibParams
    opts.MaxIter=Options.MaxIter;                         % The maximum number of iterations
    opts.PopSize=20;                                      % The population size
    %Note: MaxEvals will be roughly MaxIter*PopSize (not exactly due to resampling)
    opts.VerboseModulo=10;                                % Display results after every 10'th iteration
    opts.TolFun=Options.Tolerance;                        % Function tolerance
    opts.TolX=Options.Tolerance;                          % Tolerance in the parameters
    opts.Plotting='on';                                   % Dislpay plotting or not
    opts.Saving='on';                                     % Saving results
    opts.SaveFileName='./Output/Calib/CMAES_run_data.mat';% Results are saved to this file
    [CalibParamsVec, fval, counteval, exitflag] = Estimate_Model_CMAES(CalibFn,CalibParams_initvalsVec,CMAES_sigma,CMAES_Insigma,opts);
    save ./Output/Calib/CMAESoutput.mat CalibParamsVec fval counteval exitflag
    time=toe;
    disp(' ')
    fprintf('End of calibration loop(CMA-ES algorithm). Time used: %.4f \n', time)
    
    % 7.3.2 Store calibrated parameters with direct parameters together
    Calib_Params=Params;
    nCalibParamsVec=length(CalibParamsVec);
    for iCalibParam = 1:nCalibParamsVec
        Calib_Params.(ParamNamesToEstimate{iCalibParam})=CalibParamsVec(iCalibParam); % Those are the final calibrated parameters (including the direct ones, here just replace the calibrated ones)
    end
        
    % 7.3.3 Get the Markov-chain transition probability matrix, invariant probability measures and endowment grid/vector, VF, and the stationary distribution from the estimated parameters
    disp(' ')
    disp('Retrieving value function and stationary dsitribution from estimated parameters...')
    [e_grid,Gamma,gammastar_calib,gammastarfull_calib]=Create_Exog_Shock(Calib_Params);
    V0=ones(n_a,n_z); % Initial guess for the value function
    [V, Policy]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Calib_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2);
    StationaryDistr=StationaryDist(Policy,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,Options.NCores,Options.Simperiods,Options.Burnin,Options.Tolerance,Options.MaxIter);
    
    % 7.3.4 Get the model statistics
    CalibStats_Calib=ModelStatistics(Options.NSim,StationaryDistr, Policy, Phi_aprimeMatrix, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Calib_Params);
    ExtraStats_Calib=ModelStatistics_ExtraStats(StationaryDistr, Policy, n_d, n_a, n_z, d_grid, a_grid, z_grid, Calib_Params);

    % 7.3.5 Save the model estimation target values based on the final calibrated parameters
    save ./Output/Calib/Calib_Results.mat p_eqm V Policy a_grid StationaryDistr Calib_Params CalibStats_Calib ExtraStats_Calib
    
elseif Options.MSM==0
    % 7.3.6 Save the model estimation target values based on the initial calibrated parameters
    % Note: no loop, values directly taken from initial parametrization, which is a result of a previous good calibration run.
    disp(' ')
    disp('Calibration without running calibration loop. Pre-calibration values set as default calibration result.')
    Calib_Params=Params;
    CalibStats_Calib=ModelStats_initvals;
    ExtraStats_Calib=ExtraStats_initvals;
    save ./Output/Calib/Calib_Results.mat p_eqm V Policy a_grid StationaryDistr Calib_Params CalibStats_Calib ExtraStats_Calib
end 

% 7.4 Check interest rates
disp(' ')
fprintf('So we have p_eqm (calculated from K/Y from data) as %.4f \n', p_eqm(1))  % alpha*1/CalibStats_Data(1)-delta
SSvalueParamNames(1).Names={};
SSvaluesFn_1 = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
SSvaluesFn_2 = @(lab,anext,a,z,e1,e2,e3,e4) lab*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvaluesFn={SSvaluesFn_1, SSvaluesFn_2};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDistr, Policy, SSvaluesFn, Calib_Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);
p_market=Calib_Params.alpha*(SSvalues_AggVars(1)^(Calib_Params.alpha-1))*(SSvalues_AggVars(2)^(1-Calib_Params.alpha))-Calib_Params.delta;
fprintf('While p_market (calculated from marginal product of capital) as %.4f \n', p_market)
Y=(SSvalues_AggVars(1)^Calib_Params.alpha)*(SSvalues_AggVars(2)^(1-Calib_Params.alpha));
K=SSvalues_AggVars(1);
L=SSvalues_AggVars(2);
fprintf('This is because: \n')
fprintf('Y = %.4f, K = %.4f, L = %.4f \n', Y,K,L)
fprintf('alpha = %.4f, delta = %.4f \n', Calib_Params.alpha, Calib_Params.delta)
p_eqm_check=Calib_Params.alpha*1/(K/Y)-Calib_Params.delta;
fprintf('Double check of p_market (based on K/Y): %.4f \n', p_eqm_check)
p_eqm_check=Calib_Params.alpha*1/(CalibStats_Calib(1))-Calib_Params.delta;
fprintf('Double check of p_market (based on Model Statistics): %.4f \n', p_eqm_check)
fprintf('Note: Both the checks of p_market should give the same value as p_market \n')

disp(' ')
fprintf('End of calibration. \n')

%% 8. Post-calibration (General Equilibrium calculation)

disp(' ')
fprintf('Initializing post-calibration (General Equilibrium calculation)... \n')
GEPriceParamNames={'r','kappa'};

% 8.1 Calculate GE of the baseline model associated with the calibrated values.
if Options.SkipGE==0 % If on the server - solve for the GE
    % Note: The algorithm used for calculating the general equilibrium is not
    % what you would want to use normally for solving this model. It finds the
    % general equilibrium using a discrete grid on interest rates, rather than
    % just solving the fixed-point problem on interest rates directly by using
    % optimization. This is done for robustness reasons.
    
    GECounter=0; % Create and save a counter, then inside GE calculation loop one can load, increment, and resave, thus keeping count.
    save ./Output/GE/GECounter.mat GECounter
    
    % 8.1.1 Create descriptions of steady state (SS) values as functions of d_grid, a_grid, z_grid & Gamma (used to calculate the integral across the SS dist fn of the functions here defined)
    [e_grid,Gamma,gammastar,gammastarfull]=Create_Exog_Shock(Calib_Params); % Markov-chain transition probability matrix, invariant probability measures and endowment grid/vector    SSvalueParamNames(1).Names={};
    SSvaluesFn_1 = @(lab,anext,a,z) a; %K
    SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
    SSvaluesFn_2 = @(lab,anext,a,z,e1,e2,e3,e4) lab*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
    SSvalueParamNames(3).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
    SSvaluesFn_TaxRevenue = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa);
    SSvalueParamNames(4).Names={'J','omega'};
    SSvaluesFn_Pensions = @(lab,anext,a,z,J,omega) omega*(z>J); % If the agent is retired, he/she earns pension omega (otherwise it is zero).
    SSvaluesFn={SSvaluesFn_1,SSvaluesFn_2,SSvaluesFn_TaxRevenue,SSvaluesFn_Pensions};
    
    % 8.1.2 Define functions for the General Equilibrium conditions
    % The GE equations typically include Market Clearance condtions, but often also things such as Government Budget balance.
    % Note: the closer to zero the value given by the function is, the closer the GE condition is to being met.

        % (1) Requirement that the interest rate equals the marginal product of capital: i.e. r = alpha*K^(alpha-1)*L^(1-alpha)-delta
    GeneralEqmEqnParamNames(1).Names={'alpha','delta'}; % Names of parameters that are needed to evaluate the GeneralEqmEqns (parameters not determined as part of GE)
    GeneralEqmEqn_1 = @(AggVars,p,alpha,delta) p(1) - (alpha*(AggVars(1)^(alpha-1))*(AggVars(2)^(1-alpha))-delta); 

        % (2) Government budget balance: i.e. G/Y + Tr/Y = T/Y or G/Y + Tr/Y - T/Y = 0
    GeneralEqmEqnParamNames(2).Names={'G','alpha'};
    GeneralEqmEqn_2 = @(AggVars,p,G,alpha) G + AggVars(4)/((AggVars(1)^alpha)*(AggVars(2)^(1-alpha))) - AggVars(3)/((AggVars(1)^alpha)*(AggVars(2)^(1-alpha))); 
    % Note: The roles of 'kappa', which is contained in p(2), is already captured in the total revenue of income taxes (AggVars(3)) and in total income (Y)).

    GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2};
       
    % 8.1.3 Grid of prices or solvers of the model General Equilibrium -
    % This is irrelevant for this type of GE finding algorithm!
    % Feasible values for interest rate
    r_grid=linspace(p_eqm(1)-0.002,p_eqm(1)+0.002,n_r)';
    % Feasible values for tax rate (Search close to given value for 'kappa')
    kappa_grid=linspace(0.90*Calib_Params.kappa,1.10*Calib_Params.kappa,n_kappa)'; % 
    p_grid=[r_grid; kappa_grid];
    fprintf('Length of GE intrest rate solver grid: %.f \n', n_p(1))
    fprintf('Length of GE budget balance solver grid: %.f \n', n_p(2))
        
    % 8.1.4 Find the competitive equilibrium in an Heterogeneous Agents setup
    disp(' ')
    disp('Starting to solve General Equilibrium of the baseline model with VFI (Heterogeneous Agents). It will take a while...')
    tic;
    V0=ones(n_a,n_z); % Initial guess for the VF
    [p_eqm,MarketClearance]=HeteroAgentStationaryEqm(V0, n_d, n_a, n_z,Gamma, d_grid, a_grid, z_grid,Phi_aprimeMatrix, UtilityFn, SSvaluesFn, GeneralEqmEqns, Calib_Params, DiscountFactorParamNames, UtilityFnParamNames, SSvalueParamNames, PhiaprimeParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, Options.NCores, Options.Simperiods, Options.Burnin, Options.Tolerance, Options.MaxIter, Options.Howards, Options.Howards2);
   
    % 8.1.5 Evaluate a few objects at the equilibrium
    disp(' ')
    disp(' ')
    disp('Evaluating a few objects at General Equilibrium...')
    Calib_Params.r=p_eqm(1);
    Calib_Params.kappa=p_eqm(2);
    Calib_Params.w=(1-Calib_Params.alpha)*(((Calib_Params.r+Calib_Params.delta)/(Calib_Params.alpha))^(Calib_Params.alpha/(Calib_Params.alpha-1)));
    [V, Policy]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Calib_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Options.Tolerance, Options.Howards, Options.Howards2);
    StationaryDistr=StationaryDist(Policy,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,Options.NCores,Options.Simperiods,Options.Burnin,Options.Tolerance,Options.MaxIter);
        
    % 8.1.6 Get the baseline model statistics
    ModelStats_GE=ModelStatistics(Options.NSim,StationaryDistr, Policy, Phi_aprimeMatrix, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Calib_Params);
    ExtraStats_GE=ModelStatistics_ExtraStats(StationaryDistr, Policy, n_d, n_a, n_z, d_grid, a_grid, z_grid, Calib_Params);
    GEtime=toc;
    
    % 8.1.7 Save the results of the GE calculation
    save Output/GE/GE_results.mat p_eqm MarketClearance V Policy a_grid StationaryDistr Calib_Params ModelStats_GE ExtraStats_GE GEtime
    disp(' ')
    fprintf('Run time of General Equilibrium calculation algorithm was: %.4f \n', GEtime)

% 8.2 Load GE results and statistics from previous findings in order not to re-run the GE-calculation subroutine
elseif Options.SkipGE==1
    load Output/GE/GE_results.mat p_eqm MarketClearance V Policy a_grid StationaryDistr Calib_Params ModelStats_GE ExtraStats_GE GEtime
    disp(' ')
    disp('Baseline General Equilibrium computation not activated. Retrieve baseline General Equilibrium results from previous findings.')
end

% 8.3 Check interest rates
disp(' ')
fprintf('So we have p_eqm (calculated with Het. Agents calculation of GE) as %.4f \n', p_eqm(1))
SSvalueParamNames(1).Names={};
SSvaluesFn_1 = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
SSvaluesFn_2 = @(lab,anext,a,z,e1,e2,e3,e4) lab*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvaluesFn={SSvaluesFn_1, SSvaluesFn_2};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDistr, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);
p_market=Params.alpha*(SSvalues_AggVars(1)^(Params.alpha-1))*(SSvalues_AggVars(2)^(1-Params.alpha))-Params.delta;
fprintf('While p_market (calculated from marginal product of capital) as %.4f \n', p_market)
Y=(SSvalues_AggVars(1)^Params.alpha)*(SSvalues_AggVars(2)^(1-Params.alpha));
K=SSvalues_AggVars(1);
L=SSvalues_AggVars(2);
fprintf('This is because: \n')
fprintf('Y = %.4f, K = %.4f, L = %.4f \n', Y,K,L)
fprintf('alpha = %.4f, delta = %.4f \n', Params.alpha, Params.delta)
p_eqm_check=Params.alpha*1/(K/Y)-Params.delta;
fprintf('Double check of p_market (based on K/Y): %.4f \n', p_eqm_check)
p_eqm_check=Params.alpha*1/(ModelStats_GE(1))-Params.delta;
fprintf('Double check of p_market (based on Model Statistics): %.4f \n', p_eqm_check)
fprintf('Note: Both the checks of p_market should give the same value as p_market \n')

% 8.4 Retrieve the value of government expenditure 'G' and resave results of calibration
Calib_Params.G=Params.G*Y;
save Output/GE/GE_results.mat p_eqm MarketClearance V Policy a_grid StationaryDistr Calib_Params ModelStats_GE ExtraStats_GE GEtime

disp(' ')
fprintf('End of post-calibration (General Equilibrium calculation). \n')
