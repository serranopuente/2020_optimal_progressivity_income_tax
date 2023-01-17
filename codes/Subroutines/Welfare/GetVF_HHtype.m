function AggValue = GetVF_HHtype(V, StationaryDistr, Policy, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames, evaluatefn, lb_percentile, ub_percentile, shock_type)
% It returns the VF of a HH type group: by percentile group of sorting variable (welath, income, VF, etc...), and by type of shock (s=1, s=2, etc...)

Outcome =  SSvalues_Indiv_Vars_Ev(V, StationaryDistr, Policy, n_d, n_a, n_z, d_grid, a_grid, z_grid, Params, SSvaluesFn, SSvalueParamNames);
if shock_type == 0
if evaluatefn == 1
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=lb_percentile).*(Outcome(:,size(Outcome,2))<ub_percentile)]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    AggValue(1)=sum(sum((Outcome(:,4)).*Outcome(:,size(Outcome,2)))); % Compute aggregate value function for individuals in the percentile of interest
elseif evaluatefn == 0 % if sorting is done by VF
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)-1); % Sort by variable of interest
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=lb_percentile).*(Outcome(:,size(Outcome,2))<ub_percentile)]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    AggValue=sum(sum((Outcome(:,4)).*Outcome(:,size(Outcome,2)))); % Compute aggregate value function for individuals in the percentile of interest
end
else
if evaluatefn == 1
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)); % Sort by variable of interest
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=lb_percentile).*(Outcome(:,size(Outcome,2))<ub_percentile).*(Outcome(:,2)==shock_type)]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    AggValue(1)=sum(sum((Outcome(:,4)).*Outcome(:,size(Outcome,2)))); % Compute aggregate value function for individuals in the percentile of interest
elseif evaluatefn == 0 % if sorting is done by VF
    Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2)))=Outcome(:,(size(Outcome,2)-(length(SSvaluesFn)-1)):(size(Outcome,2))).*(Outcome(:,3)>0); % Discard points of the grid with no presence of agents
    Outcome=sortrows(Outcome,size(Outcome,2)-1); % Sort by variable of interest
    Outcome=[Outcome cumsum(Outcome(:,3))/sum(Outcome(:,3))]; % Get percentiles of variable of interest
    Outcome=[Outcome Outcome(:,3).*(Outcome(:,size(Outcome,2))>=lb_percentile).*(Outcome(:,size(Outcome,2))<ub_percentile).*(Outcome(:,2)==shock_type)]; % Get stationary distribution with 0 out of percentile of interest
    Outcome=[Outcome Outcome(:,size(Outcome,2))/sum(sum(Outcome(:,size(Outcome,2))))]; % Adjust stationary distribution just to the percentile of interest
    AggValue=sum(sum((Outcome(:,4)).*Outcome(:,size(Outcome,2)))); % Compute aggregate value function for individuals in the percentile of interest
end
    
end
end