function [PricePathNew]=TransitionPath(PricePathOld, PriceParamNames, ParamPath, PathParamNames, Parameters, DiscountFactorParamNames, Phi_aprimeKron_final, T, V_final, StationaryDist_init, ReturnFn, ReturnFnParamNames, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, SSvaluesFn,SSvalueParamNames, MarketPriceEqns, MarketPriceParamNames,Tolerance,weightscheme,oldpathweight)
%%

N_d=prod(n_d);
N_z=prod(n_z);
N_a=prod(n_a);
l_p=size(PricePathOld,2);

PricePathDist=Inf;
pathcounter=1;

V_final=reshape(V_final,[N_a,N_z]); 
AgentDist_initial=reshape(StationaryDist_init,[N_a*N_z,1]);
V=zeros(size(V_final));
PricePathNew=zeros(size(PricePathOld)); PricePathNew(T,:)=PricePathOld(T,:);
Policy=zeros(N_a,N_z);
Phi_aprimeKron=Phi_aprimeKron_final; % Might want to change this so that Phi_aprimeKron can change along the transition path.


% del otro trans
% V_final=reshape(V_final,[N_a,N_z]);
% SteadyStateDist_initial=reshape(SteadyStateDist_initial,[N_a*N_z,1]);
% V=zeros(size(V_final));
% PricePathNew=zeros(size(PricePathOld)); PricePathNew(T)=PricePathOld(T);
% PolicyIndexes=zeros(N_a,N_z);
% 
% Phi_aprimeKron=reshape(Phi_aprimeKron, [N_d,N_a*N_z,N_z]); % esto es distinto

beta=CreateVectorFromParams(Parameters, DiscountFactorParamNames);
IndexesForPathParamsInDiscountFactor=CreateParamVectorIndexes(DiscountFactorParamNames, PathParamNames);
IndexesForDiscountFactorInPathParams=CreateParamVectorIndexes(PathParamNames,DiscountFactorParamNames);
ReturnFnParamsVec=CreateVectorFromParams(Parameters, ReturnFnParamNames);
IndexesForPricePathInReturnFnParams=CreateParamVectorIndexes(ReturnFnParamNames, PriceParamNames);
IndexesForReturnFnParamsInPricePath=CreateParamVectorIndexes(PriceParamNames, ReturnFnParamNames);
IndexesForPathParamsInReturnFnParams=CreateParamVectorIndexes(ReturnFnParamNames, PathParamNames);
IndexesForReturnFnParamsInPathParams=CreateParamVectorIndexes(PathParamNames,ReturnFnParamNames);
SSvalueParamsVec=CreateVectorFromParams(Parameters, SSvalueParamNames);
IndexesForPricePathInSSvalueParams=CreateParamVectorIndexes(SSvalueParamNames, PriceParamNames);
IndexesForSSvalueParamsInPricePath=CreateParamVectorIndexes(PriceParamNames,SSvalueParamNames);
IndexesForPathParamsInSSvalueParams=CreateParamVectorIndexes(SSvalueParamNames, PathParamNames);
IndexesForSSvalueParamsInPathParams=CreateParamVectorIndexes(PathParamNames,SSvalueParamNames);
MarketPriceParamsVec=CreateVectorFromParams(Parameters, MarketPriceParamNames);
IndexesForPricePathInMarketPriceParams=CreateParamVectorIndexes(MarketPriceParamNames, PriceParamNames);
IndexesForMarketPriceParamsInPricePath=CreateParamVectorIndexes(PriceParamNames, MarketPriceParamNames);
IndexesForPathParamsInMarketPriceParams=CreateParamVectorIndexes(MarketPriceParamNames, PathParamNames);
IndexesForMarketPriceParamsInPathParams=CreateParamVectorIndexes(PathParamNames,MarketPriceParamNames);

%%
aaa=kron(pi_z,ones(N_d,1));

while PricePathDist>Tolerance
    PolicyIndexesPath=zeros(N_a,N_z,T-1); %Periods 1 to T-1

    %First, go from T-1 to 1 calculating the Value function and Optimal
    %policy function at each step. Since we won't need to keep the value
    %functions for anything later we just store the next period one in
    %Vnext, and the current period one to be calculated in V
    Vnext=V_final;
    for i=1:T-1 %so t=T-i
        
        if ~isnan(IndexesForPathParamsInDiscountFactor)
                beta(IndexesForPathParamsInDiscountFactor)=ParamPath(T-i,IndexesForDiscountFactorInPathParams); % This step could be moved outside all the loops
        end
        if ~isnan(IndexesForPricePathInReturnFnParams)
                ReturnFnParamsVec(IndexesForPricePathInReturnFnParams)=PricePathOld(T-i,IndexesForReturnFnParamsInPricePath);
        end
        if ~isnan(IndexesForPathParamsInReturnFnParams)
                ReturnFnParamsVec(IndexesForPathParamsInReturnFnParams)=ParamPath(T-i,IndexesForReturnFnParamsInPathParams); % This step could be moved outside all the loops by using BigReturnFnParamsVec idea
        end
        ReturnMatrix=CreateReturnFnMatrix(ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid, ReturnFnParamsVec);
        fprintf('%i out of %i \n',i,T-1)

        EV=zeros(N_d*N_z,N_z);
        for zprime_c=1:N_z
            EV(:,zprime_c)=Vnext(Phi_aprimeKron(:,:,zprime_c),zprime_c); %(d,z')
        end
        EV=EV.*aaa;
        EV(isnan(EV))=0; %multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
        EV=reshape(sum(EV,2),[N_d,1,N_z]);

        for z_c=1:N_z % Can probably eliminate this loop and replace with a matrix multiplication operation thereby making it faster
            entireRHS=ReturnMatrix(:,:,z_c)+beta*EV(:,z_c)*ones(1,N_a,1);

            %Calc the max and it's index
            [Vtemp,maxindex]=max(entireRHS,[],1);
            V(:,z_c)=Vtemp;
            Policy(:,z_c)=maxindex;
        end

        PolicyIndexesPath(:,:,T-i)=Policy;
        Vnext=V;
    end


    %Now we have the full PolicyIndexesPath, we go forward in time from 1
    %to T using the policies to update the agents distribution generating a
    %new price path
    %Call AgentDist the current periods distn and AgentDistnext
    %the next periods distn which we must calculate
    AgentDist=AgentDist_initial;
    for i=1:T-1
        %Get the current optimal policy
        Policy=PolicyIndexesPath(:,:,i);

        % optaprime is here replaced by Phi_of_Policy, which is a different shape
        Phi_of_Policy=zeros(N_a,N_z,N_z); %a'(a,z',z)
        for z_c=1:N_z
            Phi_of_Policy(:,:,z_c)=Phi_aprimeKron(Policy(:,z_c),:,z_c);
        end
        Ptemp=zeros(N_a,N_a*N_z*N_z);
        Ptemp(reshape(permute(Phi_of_Policy,[2,1,3]),[1,N_a*N_z*N_z])+N_a*(0:1:N_a*N_z*N_z-1))=1;
        Ptran=kron(pi_z',ones(N_a,N_a)).*reshape(Ptemp,[N_a*N_z,N_a*N_z]);
        AgentDistnext=Ptran*AgentDist;

        p=PricePathOld(i,:);
        
        if ~isnan(IndexesForPricePathInSSvalueParams)
            SSvalueParamsVec(IndexesForPricePathInSSvalueParams)=PricePathOld(i,IndexesForSSvalueParamsInPricePath);
        end
        if ~isnan(IndexesForPathParamsInSSvalueParams)
            SSvalueParamsVec(IndexesForPathParamsInSSvalueParams)=ParamPath(i,IndexesForSSvalueParamsInPathParams); % This step could be moved outside all the loops by using BigReturnFnParamsVec idea
        end
        PolicyTemp=UnKronPolicyIndexes(Policy, n_d, n_a, n_z,0);
        SSvalues_AggVars=SSvalues_AggVars_vec(AgentDist, PolicyTemp, SSvaluesFn, SSvalueParamsVec, n_d, n_a, n_z, d_grid, a_grid, z_grid);

        if ~isnan(IndexesForPricePathInMarketPriceParams)
            MarketPriceParamsVec(IndexesForPricePathInMarketPriceParams)=PricePathOld(i,IndexesForMarketPriceParamsInPricePath);
        end
        if ~isnan(IndexesForPathParamsInMarketPriceParams)
            MarketPriceParamsVec(IndexesForPathParamsInMarketPriceParams)=ParamPath(i,IndexesForMarketPriceParamsInPathParams); % This step could be moved outside all the loops by using BigReturnFnParamsVec idea
        end
        %An easy way to get the new prices is just to call MarketClearance
        %and then adjust it for the current prices
        MarketPriceParamssCell=num2cell(MarketPriceParamsVec);
        for j=1:length(MarketPriceEqns)

            % When using negative powers matlab will often return complex
            % numbers, even if the solution is actually a real number. I
            % force converting these to real, albeit at the risk of missing problems
            % created by actual complex numbers.
            PricePathNew(i,j)=real(MarketPriceEqns{j}(SSvalues_AggVars,p, MarketPriceParamssCell{:}));

        end

        AgentDist=AgentDistnext;
    end

    %See how far apart the price paths are
    PricePathDist=max(abs(reshape(PricePathNew(1:T-1,:)-PricePathOld(1:T-1,:),[numel(PricePathOld(1:T-1,:)),1])));
    
    % PricePathDist=sum(abs(reshape(PricePathNew-PricePathOld,[numel(PricePathOld),1])));

    
    %Notice that the distance is always calculated ignoring the time t=1 &
    %t=T periods, as these needn't ever converges
    disp(' ')
    disp('    Old        Old       New        New')
    [PricePathOld,PricePathNew]
   

    %Set price path to be 9/10ths the old path and 1/10th the new path (but
    %making sure to leave prices in periods 1 & T unchanged).
    if weightscheme==1 % Just a constant weighting
        PricePathOld(1:T-1,:)=oldpathweight*PricePathOld(1:T-1)+(1-oldpathweight)*PricePathNew(1:T-1,:);
    elseif weightscheme==2 % A exponentially decreasing weighting on new path from (1-oldpathweight) in first period, down to 0.1*(1-oldpathweight) in T-1 period.
        % I should precalculate these weighting vectors
        PricePathOld(1:T-1,:)=((oldpathweight+(1-exp(linspace(0,log(0.2),T-1)))*(1-oldpathweightt))'*ones(1,l_p)).*PricePathOld(1:T-1,:)+((exp(linspace(0,log(0.2),T-1)).*(1-oldpathweight))'*ones(1,l_p)).*PricePathNew(1:T-1,:);
    elseif weightscheme==3
        if (pathcounter*3)<T-1
            PricePathOld(1:(pathcounter*3),:)=oldpathweight*PricePathOld(1:(pathcounter*3),:)+(1-oldpathweight)*PricePathNew(1:(pathcounter*3),:);
        else
            PricePathOld(1:T-1,:)=oldpathweight*PricePathOld(1:T-1,:)+(1-oldpathweight)*PricePathNew(1:T-1,:);
        end
    elseif weightscheme==4 % Combines weightscheme 2 & 3
        if (pathcounter*3)<T-1
            PricePathOld(1:(pathcounter*3),:)=((oldpathweight+(1-exp(linspace(0,log(0.2),pathcounter*3)))*(1-oldpathweight))'*ones(1,l_p)).*PricePathOld(1:(pathcounter*3),:)+((exp(linspace(0,log(0.2),pathcounter*3)).*(1-oldpathweight))'*ones(1,l_p)).*PricePathNew(1:(pathcounter*3),:);
        else
            PricePathOld(1:T-1,:)=((oldpathweight+(1-exp(linspace(0,log(0.2),T-1)))*(1-oldpathweight))'*ones(1,l_p)).*PricePathOld(1:T-1,:)+((exp(linspace(0,log(0.2),T-1)).*(1-oldpathweight))'*ones(1,l_p)).*PricePathNew(1:T-1,:);
        end
    end
    
    % PricePathOld(2:T-1)=0.9*PricePathOld(2:T-1)+0.1*PricePathNew(2:T-1);
    
    
    
    TransPathConvergence=PricePathDist/Tolerance; %So when this gets to 1 we have convergence (uncomment when you want to see how the convergence isgoing)
   
    disp(' ')
    fprintf('Number of iterations on transition path: %i \n',pathcounter)
    fprintf('Current distance to convergence: %.2f (convergence when reaches 1) \n',TransPathConvergence) %So when this gets to 1 we have convergence (uncomment when you want to see how the convergence isgoing)
    
    save ./Output/Transition/TransPathConv.mat TransPathConvergence pathcounter

    %     if pathcounter==1
    %         save ./SavedOutput/FirstTransPath.mat V_final V PolicyIndexesPath PricePathOld PricePathNew
    %     end

    pathcounter=pathcounter+1;
end

end