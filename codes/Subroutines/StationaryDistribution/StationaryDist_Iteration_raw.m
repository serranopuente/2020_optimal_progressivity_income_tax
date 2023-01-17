function StationaryDistKron=StationaryDist_Iteration_raw(StationaryDistKron,PolicyKron,Phi_aprimeKron,N_a,N_z,Gamma,Tolerance,MaxIter)
% First, generate the transition matrix P=g of Q (the convolution of the 
% optimal policy function and the transition fn for exogenous shocks)

Phi_of_Policy=zeros(N_a,N_z,N_z); % Phi_of_Policy is a'(a,z,z')
for z_c=1:N_z % z
    Phi_of_Policy(:,z_c,:)=Phi_aprimeKron(PolicyKron(:,z_c),z_c,:);
end
Ptemp=zeros(N_a,N_a*N_z*N_z); % (a,z)-by-(a',z')
Ptemp(reshape(permute(Phi_of_Policy,[3,1,2]),[1,N_a*N_z*N_z])+N_a*((0:1:N_a*N_z*N_z-1)))=1;
Ptran=kron(Gamma',ones(N_a,N_a)).*reshape(Ptemp,[N_a*N_z,N_a*N_z]);

SteadyStateDistKronOld=zeros(N_a*N_z,1);
SScurrdist=sum(abs(StationaryDistKron-SteadyStateDistKronOld));
SScounter=0;

fprintf('Debugging Stationary Distr. (Iteration) - SS current distance of the distribution: %.2f (should be 1 prior to starting) \n',SScurrdist )

while SScurrdist>Tolerance && (100*SScounter)<50*MaxIter

    for jj=1:100
        StationaryDistKron=Ptran*StationaryDistKron; % No point checking distance every single iteration. Do 100, then check.
    end

    SteadyStateDistKronOld=StationaryDistKron;
    StationaryDistKron=Ptran*StationaryDistKron;
    SScurrdist=sum(abs(StationaryDistKron-SteadyStateDistKronOld));

    SScounter=SScounter+1;
    if rem(SScounter,50)==0
        SScounter
        SScurrdist
    end
end
    
    
if ~(SScounter<5*10^4)
    disp('WARNING: Steady State stopped due to reaching 50*MaxIter, this might be causing a problem')
end 

end