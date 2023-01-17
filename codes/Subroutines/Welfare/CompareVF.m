function AggVFDist = CompareVF(V0, CEV, V_ref, StationaryDistr_ref, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2, NCores, Simperiods, Burnin, MaxIter, comparison_type)

% It evaluates CEV (Consumption Equivalent Variation) between baseline
% economy and reformed economy in order to get the welfare gains / losses
% of the reform.

% Evaluate VF and Stationary Distr. with CEV
[V_bench, Policy_bench]=ValueFnIter_CEV(V0, CEV, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprimeMatrix, UtilityFn, Params, DiscountFactorParamNames, UtilityFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2);
StationaryDistr_bench=StationaryDist(Policy_bench,Phi_aprimeMatrix,n_d,n_a,n_z,Gamma,NCores,Simperiods,Burnin,Tolerance,MaxIter);

% Compare Agg. Welfare of reformed economy and baseline economy with CEV
Agg_Welfare_bench = sum(sum(V_bench.*StationaryDistr_bench));
Agg_Welfare_ref = sum(sum(V_ref.*StationaryDistr_ref));

% Compute the distance (this should be minimized)
AggVFDist = abs(Agg_Welfare_bench - Agg_Welfare_ref);

% Display status of the loop
disp(' ')
disp('Currently:')
load ./Output/Reform/ReformCounter.mat ReformCounter MaxReforms
load ./Output/Reform/IterCEVCounter.mat IterCEVCounter
IterCEVCounter=IterCEVCounter+1;
save ./Output/Reform/IterCEVCounter.mat IterCEVCounter
if comparison_type == 0
    fprintf('Calculating CEV to determine welfare gains of the reformed economy %.f out of %.f, iteration %.f... \n', ReformCounter, MaxReforms, IterCEVCounter);
    fprintf('For CEV = %.4f, distance between aggregate welfare of the baseline economy and reformed economy is: %.7f \n', CEV, AggVFDist);
elseif comparison_type == 1
    fprintf('Calculating CEV to determine welfare gains of a hypothetical economy ignoring changes in the distribution of HH after reform %.f out of %.f, iteration %.f... \n', ReformCounter, MaxReforms, IterCEVCounter);
    fprintf('For CEV = %.4f, distance between aggregate welfare of the baseline economy and hypothetical economy is: %.7f \n', CEV, AggVFDist);
elseif comparison_type == 2
    fprintf('Calculating CEV to determine welfare gains of a hypothetical economy ignoring changes in the distribution of HH and changes in size of the economy/equilibrium prices after reform %.f out of %.f, iteration %.f... \n', ReformCounter, MaxReforms, IterCEVCounter);
    fprintf('For CEV = %.4f, distance between aggregate welfare of the baseline economy and hypothetical economy is: %.7f \n', CEV, AggVFDist);
end    
    
end