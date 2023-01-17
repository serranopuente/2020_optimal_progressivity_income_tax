function [V, Policy]=ValueFnIter(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, Gamma, Phi_aprime, ReturnFn, Parameters, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, Tolerance, Howards, Howards2)

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

%% Check if the Markov transition probability matrix is so (positive values and rows summing up to one)
if     min(min(Gamma))<0
       fprintf('ERROR: Problem with Gamma in ValueFnIter: min(min(Gamma))<0 \n')
       min(min(Gamma))
       dbstack
       return
elseif max(abs((sum(Gamma,2))-1)) > 10^(-13)
       fprintf('WARNING: Problem with Gamma in ValueFnIter: rows do not sum up to one \n')
       max(sum(Gamma,2))
       min(sum(Gamma,2))
       dbstack
       return
end

% Create a vector containing all the function parameters (in order)
ReturnFnParamsVec=CreateVectorFromParams(Parameters, ReturnFnParamNames);
DiscountFactorParamsVec=CreateVectorFromParams(Parameters, DiscountFactorParamNames);
PhiaprimeParamsVec=CreateVectorFromParams(Parameters, PhiaprimeParamNames);
    
%% CreateReturnFnMatrix creates a matrix of dimension d-by-a-by-z.
% Note: since the utility/return function is independent of time creating it once and
% then using it every iteration is good for speed, but it does use a lot of memory.
ReturnMatrix=CreateReturnFnMatrix(ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid, ReturnFnParamsVec);
    
%% The Value Function Iteration
V0Kron=reshape(V0,[N_a,N_z]);

[VKron, Policy]=ValueFnIter_raw(V0Kron, n_d,n_a,n_z, Gamma, Parameters.beta, ReturnMatrix, Phi_aprime, Tolerance, Howards, Howards2); 
    
%% Sort out Policy

V=reshape(VKron,[n_a,n_z]);
PolicyTemp=zeros(length(n_d),N_a,N_z);
    for i=1:N_a
        for j=1:N_z
            optdindexKron=Policy(i,j);
            optD=ind2sub_homemade(n_d',optdindexKron);
            PolicyTemp(:,i,j)=[optD'];
        end
    end
Policy=reshape(PolicyTemp,[length(n_d),n_a,n_z]);


end