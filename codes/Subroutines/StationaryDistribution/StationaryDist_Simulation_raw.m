function StationaryDistKron=StationaryDist_Simulation_raw(PolicyIndexesKron,Phi_aprimeKron,N_a,N_z,Gamma,NCores,Simperiods,Burnin)
% Simulates a path based on PolicyIndexes (and Phi_aprime) of length 'periods' after a burn
% in of length 'burnin' (burn-in are the initial run of points that are then
% dropped)

%%

eachsimperiods=ceil(Simperiods/NCores);
StationaryDistKron=zeros(N_a*N_z,NCores);
cumsum_Gamma=cumsum(Gamma,2);

parfor seed_c=1:NCores
        SteadyStateDistKron_seed_c=zeros(N_a*N_z,1);
        currstate=[ceil(N_a/2),ceil(N_z/2)]; %Pick a random start point on the (vectorized) (a,z) grid
        for i=1:Burnin
            optd=PolicyIndexesKron(currstate(1),currstate(2));
            [~,newcurrstate2]=max(cumsum_Gamma(currstate(2),:)>rand(1,1));
            currstate(1)=Phi_aprimeKron(optd,currstate(2),newcurrstate2);
            currstate(2)=newcurrstate2;
        end
        for i=1:eachsimperiods
            SteadyStateDistKron_seed_c(currstate(1)+(currstate(2)-1)*N_a)=SteadyStateDistKron_seed_c(currstate(1)+(currstate(2)-1)*N_a)+1;

            optd=PolicyIndexesKron(currstate(1),currstate(2));
            [~,newcurrstate2]=max(cumsum_Gamma(currstate(2),:)>rand(1,1));
            currstate(1)=Phi_aprimeKron(optd,currstate(2),newcurrstate2);
            currstate(2)=newcurrstate2;
        end
        StationaryDistKron(:,seed_c)=SteadyStateDistKron_seed_c;
end
    
StationaryDistKron=sum(StationaryDistKron,2)./(eachsimperiods*NCores);

end