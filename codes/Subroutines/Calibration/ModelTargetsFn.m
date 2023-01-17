function CalibStats_Model=ModelTargetsFn(CalibParamsVec,Params,n_d,n_a,n_z,a_grid,ReturnFn,ReturnFnParamNames,DiscountFactorParamNames,PhiaprimeParamNames,ParamNamesToEstimate,NumOfCalibStats,CalibStatNames,CalibStats_Data,NCores,NSim,Howards,Howards2,Simperiods,Burnin,Tolerance,MaxIter)

%% Take parameters from vector into structure
nCalibParamsVec=length(CalibParamsVec);

for iCalibParam = 1:nCalibParamsVec
    Params.(ParamNamesToEstimate{iCalibParam})=CalibParamsVec(iCalibParam);
end

%% Define the wage from directly identified parameters and the given interest rate
Params.w=(1-Params.alpha)*(((Params.r+Params.delta)/(Params.alpha))^(Params.alpha/(Params.alpha-1)));

%% %% Create all the other objects needed and which depend on the parameters to be calibrated
[e_grid,Gamma,Gamma_ee,gammastar,gammastarfull]=Create_Exog_Shock(Params); % Transition probability matrix
l_grid=linspace(0,Params.elle,n_d(1))'; % Grid on endogeneous state - labor choice
d_grid=[l_grid; a_grid]; % Bring model into the notational conventions used by the algorithm
z_grid=linspace(1,2*Params.J,2*Params.J)'; %(age (& determines retirement))
Phi_aprimeMatrix=PhiaprimeMatrix(n_d,n_z); % Matrix of next period's assets
V0=ones(n_a,n_z); % Guess on the VF

%% Calculate steady state
tic;
disp(' ')
fprintf('Model Targets Fn: Solving VFI... \n')
[V, Policy]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2);

disp(' ')
fprintf('Model Targets Fn: Solving stationary distribution... \n')
StationaryDistKron=StationaryDist(Policy,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,NCores, Simperiods, Burnin, Tolerance, MaxIter);

temp=min(min(StationaryDistKron));
if temp<0
    fprintf('Model Targets Fn:: Value of min(min(StationaryDist)) =%8.2f \n', temp )
    fprintf('Model Targets Fn: Seems that there is some kind of rounding error that sometimes gives min(min(StationaryDist))<0, so have added line to forcibly override this \n')
    StationaryDistKron(StationaryDistKron<0)=0;
    StationaryDistKron=StationaryDistKron./sum(sum(StationaryDistKron));
    fprintf('Model Targets Fn: (Corrected) Value of min(min(StationaryDist)) =%8.2f \n', min(min(StationaryDistKron)) )
end
fprintf('Model Targets Fn: Total mass of stationary dist=%8.2f \n', sum(sum(StationaryDistKron)))
disp(' ')
fprintf('Retrieving model statistics... \n')

%% Model target statistics
CalibStats_Model=ModelStatistics(NSim,StationaryDistKron, Policy, Phi_aprimeMatrix, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Params);
time=toc;
disp(' ')
fprintf('Run time of calibration step: %8.4f \n',time)


%% Print out the current state of the calibration process
disp(' ')
disp('Calibration status: Current params being tested:')
CalibParamsVec=CreateVectorFromParams(Params,ParamNamesToEstimate);
CalibParamsStr=CreateParamsStrucFromParamsVec(ParamNamesToEstimate, CalibParamsVec)

C=cell(NumOfCalibStats,3);
for i=1:NumOfCalibStats
    C(i,1)=CalibStatNames(i);
    C(i,2)={CalibStats_Data(i)};
    C(i,3)={CalibStats_Model(i)};
end
disp(' ')
disp('            Columns                                           Data              Model')
C(:,:)

load ./Output/Calib/CalibCounter.mat CalibCounter
CalibCounter=CalibCounter+1;
save ./Output/Calib/CalibCounter.mat CalibCounter
disp(' ')
disp(['Current calibration count is: ', num2str(CalibCounter)])

save ./Output/Calib/Calib.mat CalibParamsStr V Policy StationaryDistKron


end