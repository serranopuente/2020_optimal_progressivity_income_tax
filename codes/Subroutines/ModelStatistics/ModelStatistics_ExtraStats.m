function ExtraModelStatistics=ModelStatistics_ExtraStats(StationaryDist, Policy, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params)
% Returns statistics calculated from the steady state of the model, which are not needed as part of the calibration.
ExtraModelStatistics=zeros(24,1);

I_retired=[zeros(Params.J,1);ones(Params.J,1)];

%% Many of our statistics are based on ratios of aggregate variables. We start therefore by calculating the values of a number of aggregate variables that we will need.
SSvalueParamNames(1).Names={};
SSvaluesFn_K = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'delta'};
SSvaluesFn_I = @(lab,anext,a,z,delta) anext-a*(1-delta); %I
SSvalueParamNames(3).Names={'e1','e2','e3','e4'};
SSvaluesFn_L = @(lab,anext,a,z,e1,e2,e3,e4) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvalueParamNames(4).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'}; % SSvaluesFn_Consumption
SSvaluesFn_Consumption = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) ConsumptionFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa); % C
SSvalueParamNames(5).Names={};
SSvalue_K_retired = @(lab,anext,a,z) a*I_retired(z); %Assets held by the retired
SSvaluesFn={SSvaluesFn_K,SSvaluesFn_I,SSvaluesFn_L,SSvaluesFn_Consumption,SSvalue_K_retired};
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDist, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);
Y=(SSvalues_AggVars(1)^Params.alpha)*(SSvalues_AggVars(3)^(1-Params.alpha)); % Y (Production technology function)

%% Some model statistics from aggregates without further calculation
ExtraModelStatistics(1)=SSvalues_AggVars(2)/Y; % I/Y
ExtraModelStatistics(2)=SSvalues_AggVars(4); % C

%% Distributional of consumption and earnings
SSvalueParamNames(1).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'}; % SSvaluesFn_Consumption
SSvaluesFn_Consumption = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) ConsumptionFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa); % C
SSvalueParamNames(2).Names={'e1','e2','e3','e4','w','r','J'};
SSvaluesFn_earnings = @(lab,anext,a,z,e1,e2,e3,e4,w,r,J,omega) lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r; % Earnings
SSvaluesFn={SSvaluesFn_Consumption,SSvaluesFn_earnings};
SSvalues_LorenzCurves=SSvalues_LorenzCurve(StationaryDist, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,100);
%  Gini for Consumption
ExtraModelStatistics(3)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(1,:));
%  Consumption Lorenz Curve - Quintiles (%) 
ExtraModelStatistics(4:7)=100*(SSvalues_LorenzCurves(1,[40,60,80,100])-SSvalues_LorenzCurves(1,[1,41,61,81]));
%  Consumption Lorenz Curve - 90-95, 95-99, and 99-100 (%)
ExtraModelStatistics(8:10)=100*(SSvalues_LorenzCurves(1,[95,99,100])-SSvalues_LorenzCurves(1,[90,95,99]));
%  Gini for Earnings
ExtraModelStatistics(11)=Gini_from_LorenzCurve(SSvalues_LorenzCurves(2,:));
%  Earnings Lorenz Curve - Quintiles (%) 
ExtraModelStatistics(12:15)=100*(SSvalues_LorenzCurves(2,[40,60,80,100])-SSvalues_LorenzCurves(2,[1,41,61,81]));
%  Earnings Lorenz Curve - 90-95, 95-99, and 99-100 (%)
ExtraModelStatistics(16:18)=100*(SSvalues_LorenzCurves(2,[95,99,100])-SSvalues_LorenzCurves(2,[90,95,99]));

%% Fraction of total assets held by retirees
ExtraModelStatistics(19)=SSvalues_AggVars(5)/SSvalues_AggVars(1);

%% Descriptive statistics of hours and consumption
SSvalueParamNames2(1).Names={};
SSvaluesFn_H = @(lab,anext,a,z) lab; %H
SSvalueParamNames2(2).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'}; % SSvaluesFn_Consumption
SSvaluesFn_Consumption = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) ConsumptionFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa); % C
SSvaluesFn2={SSvaluesFn_H,SSvaluesFn_Consumption};
SSvalues_MeanMedianStdDev_H_C=SSvalues_MeanMedianStdDev(StationaryDist, Policy, SSvaluesFn2, Params,SSvalueParamNames2, n_d, n_a, n_z, d_grid, a_grid, z_grid);

%% Mean, Median & Std Dev of Hours Worked
ExtraModelStatistics(20)=SSvalues_MeanMedianStdDev_H_C(1,1);
ExtraModelStatistics(21)=SSvalues_MeanMedianStdDev_H_C(1,2);
ExtraModelStatistics(22)=SSvalues_MeanMedianStdDev_H_C(1,3);

%% Ratio of cross-section coeff.s of variation of consumption and hours worked
cv_of_c=SSvalues_MeanMedianStdDev_H_C(1,3)/SSvalues_MeanMedianStdDev_H_C(1,1);
cv_of_h=SSvalues_MeanMedianStdDev_H_C(2,3)/SSvalues_MeanMedianStdDev_H_C(2,1);
ExtraModelStatistics(23)=cv_of_c/cv_of_h;

%% Average effective tax rate
SSvalueParamNames3(1).Names={'J','r','w','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
SSvaluesFn_effectivetax = @(lab,anext,a,z,J,r,w,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) (TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa))/(lab.*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4))*w+a*r+omega.*(z>J));
SSvaluesFn3={SSvaluesFn_effectivetax};
SSvalues_MeanMedianStdDev_effectivetax=SSvalues_MeanMedianStdDev(StationaryDist, Policy, SSvaluesFn3, Params,SSvalueParamNames3, n_d, n_a, n_z, d_grid, a_grid, z_grid);
ExtraModelStatistics(24)=SSvalues_MeanMedianStdDev_effectivetax(1,1);

end