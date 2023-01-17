function StationaryDist=StationaryDist(Policy,Phi_aprimeKron,n_d,n_a,n_z,Gamma,NCores,Simperiods,Burnin,Tolerance,MaxIter)
N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

%%

tempPolicy=reshape(Policy,[length(n_d),N_a,N_z]); % First dim indexes the optimal choice for d and aprime rest of dimensions a,z
PolicyKron=zeros(N_a,N_z);
for i1=1:N_a
    for i2=1:N_z
        PolicyKron(i1,i2)=sub2ind_homemade([n_d],tempPolicy(:,i1,i2));
    end
end

StationaryDistKron=StationaryDist_Simulation_raw(PolicyKron,Phi_aprimeKron,N_a,N_z,Gamma,NCores,Simperiods,Burnin);
fprintf('Debugging Stationary Distr. (before iterate) - total distribution mass: %.2f \n', sum(sum(StationaryDistKron)))

StationaryDistKron=StationaryDist_Iteration_raw(StationaryDistKron,PolicyKron,Phi_aprimeKron,N_a,N_z,Gamma,Tolerance,MaxIter);
fprintf('Debugging Stationary Distr. (after iterate) - total distribution mass: %.2f \n', sum(sum(StationaryDistKron)))

StationaryDist=reshape(StationaryDistKron,[n_a,n_z]);

end