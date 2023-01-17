function ModelStatistics=ModelStatistics(NSim,StationaryDist, Policy, Phi_aprimeKron, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Params)
% Returns 22 statistics calculated from the steady state of the model
% Almost all of the runtime is in calculating the ratio of earnings old to young and intergenerational earnings (especially the first of these).

ModelStatistics=zeros(22,1);

%% Many of our statistics are based on ratios of aggregate variables. We start therefore by calculating the values of a number of aggregate variables that we will need.

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
SSvaluesFn={SSvaluesFn_K,SSvaluesFn_L,SSvaluesFn_H,SSvaluesFn_TaxRevenue,SSvaluesFn_Pensions};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDist, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);
Y=(SSvalues_AggVars(1)^Params.alpha)*(SSvalues_AggVars(2)^(1-Params.alpha)); % Y (Production technology function)
r=Params.alpha*(SSvalues_AggVars(1)^(Params.alpha-1))*(SSvalues_AggVars(2)^(1-Params.alpha))-Params.delta;

disp(' ');
disp('Model aggregate statistics:');
fprintf('Y = %.4f \n', Y)
fprintf('r in the market = %.4f \n', r)
fprintf('alpha = %.4f \n', Params.alpha)
fprintf('delta = %.4f \n', Params.delta)
fprintf('K = %.4f \n', SSvalues_AggVars(1))
fprintf('L = %.4f \n', SSvalues_AggVars(2))

%% Some model statistics from aggregates without further calculation
ModelStatistics(1)=SSvalues_AggVars(1)/Y; % K/Y - Capital-output ratio
ModelStatistics(2)=100*(SSvalues_AggVars(4)-SSvalues_AggVars(5))/Y; % G/Y:  Government-expenditure to output ratio, G=T-Tr=Tax Revenue - Pensions
ModelStatistics(3)=100*SSvalues_AggVars(5)/Y; % Tr/Y Transfers To OutputRatio
ModelStatistics(4)=100*SSvalues_AggVars(3)/Params.elle; % h Share of disposable time allocated to market

%% Ratio of earnings of 40 year olds to 20 year olds.
e=[Params.e1,Params.e2,Params.e3,Params.e4,0,0,0,0];
ModelStatistics(5)=RatioEarningsOldYoung(NSim, StationaryDist, Policy, Phi_aprimeKron, n_d,n_a,n_z, d_grid, Gamma, e,Params.w,Params.J); % Ratio of earnings of 40 year olds to 20 year olds
%% Intergenerational correlation coefficient. This is quite complicated to calculate and so required a dedicated script.
ModelStatistics(6)=IntergenerationalEarnings(NSim,StationaryDist, Policy, Phi_aprimeKron, n_d,n_a,n_z,d_grid, Gamma, e,Params.w,Params.J); % Cross sectional correlation of incomes between fathers and sons

%% Distributional statistics of income and wealth
SSvalueParamNames(1).Names={'e1','e2','e3','e4','w','r','J','omega'}; % income
SSvaluesFn_income = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J);
SSvalueParamNames(2).Names={}; % K
SSvaluesFn={SSvaluesFn_income,SSvaluesFn_K};
SSvalues_LorenzCurves=SSvalues_LorenzCurve(StationaryDist, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,100);
%  Gini for Income
ModelStatistics(7)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(1,:));
%  Income Lorenz Curve - Quintiles (%) 
ModelStatistics(9:12)=100*(SSvalues_LorenzCurves(1,[40,60,80,100])-SSvalues_LorenzCurves(1,[1,41,61,81]));
%  Income Lorenz Curve - 90-95, 95-99, and 99-100 (%)
ModelStatistics(17:19)=100*(SSvalues_LorenzCurves(1,[95,99,100])-SSvalues_LorenzCurves(1,[90,95,99]));
%  Gini for Wealth
ModelStatistics(8)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(2,:));
%  Wealth Lorenz Curve - Quintiles (%)
ModelStatistics(13:16)=100*(SSvalues_LorenzCurves(2,[40,60,80,100])-SSvalues_LorenzCurves(2,[1,41,61,81]));
%  Wealth Lorenz Curve - 90-95, 95-99, and 99-100 (%)
ModelStatistics(20:22)=100*(SSvalues_LorenzCurves(2,[95,99,100])-SSvalues_LorenzCurves(2,[90,95,99]));
% Adjust income and wealth shares to make them to sum up to 1 - this is due
% to calculus error that mak it impossible to sum up to 1 directly (it
% normally results in a sum of 98 instead of 100)
TotIncomeShares=sum(ModelStatistics(9:12));
ModelStatistics(9:12)=(ModelStatistics(9:12))/TotIncomeShares*100;
TotWealthShares=sum(ModelStatistics(13:16));
ModelStatistics(13:16)=(ModelStatistics(13:16))/TotWealthShares*100;

%% To stop the code from going out of bounds in the sense of negative values of Gamma
if min(Gamma>=0)==0
    ModelStatistics=-Inf*ones(22,1);
    disp('Bad calibration: Gamma matrix with some negative elements')
    Gamma
    sum(Gamma,2)
end

end