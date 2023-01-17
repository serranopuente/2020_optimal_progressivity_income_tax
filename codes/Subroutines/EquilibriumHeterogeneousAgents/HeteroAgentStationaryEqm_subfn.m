function GeneralEqmConditions=HeteroAgentStationaryEqm_subfn(p, V0Kron, n_d, n_a, n_z, Gamma, d_grid, a_grid, z_grid, Phi_aprimeKron, ReturnFn, SSvaluesFn, GeneralEqmEqns, Parameters, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, PhiaprimeParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, NCores, Simperiods, Burnin, Tolerance, MaxIter, Howards, Howards2)

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

%% 
for j=1:length(GEPriceParamNames)
    Parameters.(GEPriceParamNames{j})=p(j);
end

[~, Policy]=ValueFnIter(V0Kron, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeKron, ReturnFn, Parameters, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2);

% Calculate the Steady-state distance (given this price) and use it to assess market clearance
StationaryDistKron=StationaryDist(Policy,Phi_aprimeKron,n_d,n_a,n_z,Gamma,NCores,Simperiods,Burnin,Tolerance,MaxIter);
SSvalues_AggVars=SSvalues_AggVars_Ev(StationaryDistKron, Policy, SSvaluesFn, Parameters, SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);

% use of real() is a hack that could disguise errors, but I couldn't
% find why matlab was treating output as complex
GeneralEqmConditionsVec=real(GeneralEqmConditions_Ev(SSvalues_AggVars,p, GeneralEqmEqns, Parameters,GeneralEqmEqnParamNames));

GeneralEqmConditions=sum(GeneralEqmConditionsVec.^2); 
GeneralEqmConditions=gather(GeneralEqmConditions);

disp(' ');
fprintf('Current Aggregates: \n')
fprintf('K = %.4f \n', SSvalues_AggVars(1))
fprintf('L = %.4f \n', SSvalues_AggVars(2))
fprintf('Income Tax Revenue = %.4f \n', SSvalues_AggVars(3))
fprintf('Pensions = %.4f \n', SSvalues_AggVars(4))
disp(' ');
fprintf('Current GE prices: r = %.4f, kappa = %.4f \n', p(1), p(2))
fprintf('General-eqm.-conditions-minimizing residual: \n')
fprintf('(1) Interest rate equals the marginal product of capital: %.7f \n', GeneralEqmConditionsVec(1))
fprintf('(2) Government budget balance: %.7f \n', GeneralEqmConditionsVec(2))

load ./Output/GE/GECounter.mat GECounter
GECounter=GECounter+1;
save ./Output/GE/GECounter.mat GECounter
disp(' ')
disp(['Current GE calculation count is: ', num2str(GECounter)])

end