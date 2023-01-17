function AggVFDist = CompareVF_Indiv(V0, CEV, V_ref, StationaryDistr_ref, Policy_ref, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Params, Ref_Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, SSvaluesFn, SSvalueParamNames, evaluatefn, lb_percentile, ub_percentile, shock_type, Tolerance, Howards, Howards2, NCores, Simperiods, Burnin, MaxIter)

% It evaluates CEV (Consumption Equivalent Variation) between baseline
% economy and reformed economy in order to get the welfare gains / losses
% of the reform for an specific HH type - by percentile group of sorting 
% variable (welath, income, VF, etc...), and by type of shock (s=1, s=2, etc...)

% Evaluate VF and Stationary Distr. with CEV
[V_bench, Policy_bench]=ValueFnIter_CEV(V0, CEV, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2);
StationaryDistr_bench=StationaryDist(Policy_bench,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,NCores,Simperiods,Burnin,Tolerance,MaxIter);

% Compare Agg. Welfare of HH type group in reformed economy and baseline economy with CEV
AggValue_bench = GetVF_HHtype(V_bench, StationaryDistr_bench, Policy_bench, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames, evaluatefn, (lb_percentile/100), (ub_percentile/100), shock_type);
AggValue_ref = GetVF_HHtype(V_ref, StationaryDistr_ref, Policy_ref, n_d, n_a, n_z, d_grid, a_grid, z_grid, Ref_Params, SSvaluesFn, SSvalueParamNames, evaluatefn, (lb_percentile/100), (ub_percentile/100), shock_type);


% Compute the distance (this should be minimized)
AggVFDist = abs(AggValue_bench - AggValue_ref);

% Display status of the loop
disp(' ')
disp('Currently:')
load ./Output/Reform/IterCEVCounter.mat IterCEVCounter
IterCEVCounter=IterCEVCounter+1;
save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
if shock_type == 0
    fprintf('Calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f for every type of income shock, iteration %.f... \n', lb_percentile, ub_percentile, IterCEVCounter);
    fprintf('For CEV = %.4f, distance between aggregate welfare of the baseline economy and reformed economy is: %.7f \n', CEV, AggVFDist);
else
    fprintf('Calculating CEV to determine welfare gains of the optimal economy for HHs between percentiles p%.f-p%.f for income shock s = %.f, iteration %.f... \n', lb_percentile, ub_percentile, shock_type, IterCEVCounter);
    fprintf('For CEV = %.4f, distance between aggregate welfare of the baseline economy and reformed economy is: %.7f \n', CEV, AggVFDist);
end
    
end