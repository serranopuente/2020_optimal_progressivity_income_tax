function [p_eqm,GeneralEqmConditions]=HeteroAgentStationaryEqm(V0,n_d, n_a, n_z, Gamma, d_grid, a_grid, z_grid,Phi_aprimeKron, ReturnFn, SSvaluesFn, GeneralEqmEqns, Parameters, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, PhiaprimeParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, Ncores, Simperiods, Burnin, Tolerance, MaxIter, Howards, Howards2)
% It uses fminsearch to find the GE (find price vector that corresponds to MarketClearance=0)
N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

%%
zerosinphi_aprimekron=sum(sum(sum(sum(Phi_aprimeKron==0))));
disp(' ');
fprintf('If this number is not zero there is an problem with Phi_aprimeKron: %.0f', zerosinphi_aprimekron)
disp(' ');

%%
V0Kron=reshape(V0,[N_a,N_z]);
tic;
GeneralEqmConditionsFn=@(p) HeteroAgentStationaryEqm_subfn(p, V0Kron, n_d, n_a, n_z, Gamma, d_grid, a_grid, z_grid, Phi_aprimeKron, ReturnFn, SSvaluesFn, GeneralEqmEqns, Parameters, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, PhiaprimeParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, Ncores, Simperiods, Burnin, Tolerance, MaxIter, Howards, Howards2);

p0=nan(length(GEPriceParamNames),1);
for ii=1:length(GEPriceParamNames)
    p0(ii)=Parameters.(GEPriceParamNames{ii});
end

[p_eqm,GeneralEqmConditions]=fminsearch(GeneralEqmConditionsFn,p0); % Solve for the prices the GE conditions with fminsearch

HA_time=toc;
disp(' ');
fprintf('Heterogeneous Agents GE resolved. Time of calculation: %.4f', HA_time)

end