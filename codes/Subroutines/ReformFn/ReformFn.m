function ReformEqmTarget=ReformFn(p, tau, weight_BudgetBalance, Params, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, GEPriceParamNames, V0, Phi_aprimeKron, TargetsBudget, PhiaprimeParamNames, n_d, n_a, n_z, a_grid, d_grid, z_grid, Gamma, Tolerance, Howards, Howards2, NCores, Simperiods, Burnin, MaxIter)

% It gives the distance of the reformed economy to the equilibrium
% conditions. It is done to iterate on price ('r' and 'lambda') values
% until reach a new equilibrium in the reformed economy (new 'tau')

for j=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{j})=p(j);
end

Params.tau=tau; 
% VFI and stationary distribution
[~, Policy]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeKron, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2);
StationaryDistKron=StationaryDist(Policy,Phi_aprimeKron,n_d,n_a,n_z,Gamma,NCores,Simperiods,Burnin,Tolerance,MaxIter);

% Define and calculate aggregates
SSvalueParamNames(1).Names={};
SSvaluesFn_1 = @(lab,anext,a,z) a; %K
SSvalueParamNames(2).Names={'e1','e2','e3','e4'};
SSvaluesFn_2 = @(lab,anext,a,z,e1,e2,e3,e4) lab*(e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4)); % Efficiency hours worked: L
SSvalueParamNames(3).Names={'J','r','alpha','delta','omega','e1','e2','e3','e4','lambda','tau','kappa'};
SSvaluesFn_TaxRevenue = @(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa) TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa);
SSvalueParamNames(4).Names={'J','omega'};
SSvaluesFn_Pensions = @(lab,anext,a,z,J,omega) omega*(z>J); % If the agent is retired, he/she earns pension omega (otherwise it is zero).
SSvaluesFn={SSvaluesFn_1,SSvaluesFn_2,SSvaluesFn_TaxRevenue,SSvaluesFn_Pensions};
AggVars=SSvalues_AggVars_Ev(StationaryDistKron, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);

% Define target and calculate distance to it
% The target equation are the GE equations that typically include Market Clearance condtions.
% Note: the closer to zero the value given by the function is, the closer the GE condition is to being met.
    % (1) Requirement that the interest rate equals the marginal product of capital: i.e. r = alpha*K^(alpha-1)*L^(1-alpha)-delta
Target_1 = (Params.alpha*(AggVars(1)^(Params.alpha-1))*(AggVars(2)^(1-Params.alpha))-Params.delta); 

    % (2) The reform balances total tax revenues by changing average level of taxes ('lambda') since government consumption-to-output ratio is fixed
Target_2 = TargetsBudget(1)/100 +  AggVars(4)/((AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha))) - AggVars(3)/((AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha)));

    % (2) The reform balances total tax keeps Tr/Y fixed by changing normalized transfers to retirees ('omega')
Target_3 = TargetsBudget(2)/100 -  AggVars(4)/((AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha)));

ReformEqmTarget_1=abs(Target_1-p(1))/p(1);
ReformEqmTarget_2=abs(Target_2)*weight_BudgetBalance(1); % This target can be augmented, otherwise sometimes would never be achieved.
ReformEqmTarget_3=abs(Target_3)*weight_BudgetBalance(2); % This target can be augmented, otherwise sometimes would never be achieved.
ReformEqmTarget=gather(ReformEqmTarget_1+ReformEqmTarget_2++ReformEqmTarget_3);

% Display status of the loop
disp(' ')
disp('Currently:')
load ./Output/Reform/ReformCounter.mat ReformCounter MaxReforms
load ./Output/Reform/IterReformCounter.mat IterReformCounter
IterReformCounter=IterReformCounter+1;
save ./Output/Reform/IterReformCounter.mat IterReformCounter
fprintf('Reform %.f out of %.f (tau is set to %.3f), iteration %.f of the GE finding algorithm... \n', ReformCounter, MaxReforms, tau, IterReformCounter);
fprintf('r = %.4f and r in the market = %.4f \n', p(1), Target_1);
fprintf('omega = %.4f, lambda = %.4f, total tax revenue = %.2f, gov-exp-to-otuput ratio = %.2f, and transfers-to-output ratio = %.2f \n', p(3), p(2), 100*AggVars(3)/((AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha))),100*(AggVars(3)-AggVars(4))/((AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha))), 100*AggVars(4)/((AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha))));
fprintf('Reform-eqm.-target-minimizing residual: \n')
fprintf('(1) Interest rate equals the marginal product of capital: %.7f \n', ReformEqmTarget_1)
fprintf('(2) Reform balances total tax revenues keeping G/Y fixed: %.7f \n', ReformEqmTarget_2)
fprintf('(3) Reform keeps Tr/Y fixed: %.7f \n', ReformEqmTarget_3)

end