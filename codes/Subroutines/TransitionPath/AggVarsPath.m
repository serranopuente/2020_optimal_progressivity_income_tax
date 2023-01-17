function [AggVarsPath]=AggVarsPath(PricePath, n_d,n_a, n_s, pi_s, d_grid,a_grid,s_grid, beta,Phi_aprimeKron, Case2_Type,T, ParamPath, V_final, SteadyStateDist_initial, FmatrixFn, SSvaluesFn,SSvalues_AggVars_initial,SSvalues_AggVars_final)
%AggVarsPath is T+1 periods long (periods 0 (before the reforms are announced) & T are the initial and final values).

N_d=prod(n_d);
N_s=prod(n_s);
N_a=prod(n_a);

% num_a_vars=length(n_a);

V_final=reshape(V_final,[N_a,N_s]);
SteadyStateDist_initial=reshape(SteadyStateDist_initial,[N_a*N_s,1]);
V=zeros(size(V_final));
PolicyIndexes=zeros(1,N_a,N_s);

% tempPhi_aprimeKron=reshape(Phi_aprime,[num_a_vars,N_d,N_a,N_s,N_s]);
% Phi_aprimeKron=zeros([N_d,N_a,N_s,N_s]);
% for i1=1:N_d
%     for i2=1:N_a
%         for i3=1:N_s
%             for i4=1:N_s
%                 Phi_aprimeKron(i1,i2,i3,i4)=sub2ind_homemade([n_a],tempPhi_aprimeKron(:,i1,i2,i3,i4));
%             end
%         end
%     end
% end

PolicyIndexesPath=zeros(N_d,N_a,N_s,T-1); %Periods 1 to T-1

%First, go from T-1 to 1 calculating the Value function and Optimal
%policy function at each step. Since we won't need to keep the value
%functions for anything later we just store the next period one in
%Vnext, and the current period one to be calculated in V
Vnext=V_final;
for i=1:T-1 %so t=T-i
    params=ParamPath(T-i,:);
    p=PricePath(T-i,:);
    Fmatrix=reshape(FmatrixFn(p,params),[N_d,N_a,N_s]);

    for s_c=1:N_s
        %first calc the second half of the RHS (except beta)
        RHSpart2=zeros(N_d,1);
        for sprime_c=1:N_s
            if pi_s(s_c,sprime_c)~=0 %multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
                RHSpart2=RHSpart2+Vnext([Phi_aprimeKron(:,s_c,sprime_c)],sprime_c)*pi_s(s_c,sprime_c);
            end
        end
        for a_c=1:N_a
            entireRHS=Fmatrix(:,a_c,s_c)+beta*RHSpart2; %d by 1

            %then maximizing d indexes
            [V(a_c,s_c),PolicyIndexes(a_c,s_c)]=max(entireRHS,[],1);
        end
    end

    PolicyIndexesPath(:,:,:,T-i)=PolicyIndexes;
    Vnext=V;
end

%Now we have the full PolicyIndexesPath, we go forward in time from 1
%to T using the policies to generate the AggVarsPath. First though we
%put in it's initial and final values.
AggVarsPath=zeros(length(SSvaluesFn),T+1);
AggVarsPath(:,1)=SSvalues_AggVars_initial; AggVarsPath(:,T+1)=SSvalues_AggVars_final;
%Call SteadyStateDist the current periods distn and SteadyStateDistnext
%the next periods distn which we must calculate
SteadyStateDist=SteadyStateDist_initial;
for i=1:T-1
    %Get the current optimal policy
    PolicyIndexes=PolicyIndexesPath(:,:,:,i);
    %Use this to calculate the steady state distn
    P=zeros(N_a,N_s,N_a,N_s); %P(a,z,aprime,zprime)=proby of going to (a',z') given in (a,z)
    for s_c=1:N_s
        for a_c=1:N_a
            optd=PolicyIndexes(a_c,s_c);
            for sprime_c=1:N_s
                optaprime=Phi_aprimeKron(optd,s_c,sprime_c);
                P(a_c,s_c,optaprime,sprime_c)=pi_s(s_c,sprime_c)/sum(pi_s(s_c,:));
            end
        end
    end
    P=reshape(P,[N_a*N_s,N_a*N_s]);
    P=P';
    SteadyStateDistnext=P*SteadyStateDist;

    AggVarsPath(:,i)=SSvalues_AggVars_raw(SteadyStateDist, PolicyIndexes, SSvaluesFn, n_d, n_a, N_s, d_grid, a_grid,s_grid,pi_s,p); %the two zeros represent the d variables

    SteadyStateDist=SteadyStateDistnext;
end
%i=T
params=ParamPath(T,:);
p=PricePath(T,:);
Fmatrix=reshape(FmatrixFn(p,params),[N_d,N_a,N_s]);
for s_c=1:N_s
    %first calc the second half of the RHS (except beta)
    RHSpart2=zeros(N_d,1);
    for sprime_c=1:N_s
        if pi_s(s_c,sprime_c)~=0 %multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
            RHSpart2=RHSpart2+V_final([Phi_aprimeKron(:,s_c,sprime_c)],sprime_c)*pi_s(s_c,sprime_c);
        end
    end
    for a_c=1:N_a
        entireRHS=Fmatrix(:,a_c,s_c)+beta*RHSpart2; %d by 1

        %then maximizing d indexes
        [V(a_c,s_c),PolicyIndexes(a_c,s_c)]=max(entireRHS,[],1);
    end
end

AggVarsPath(:,T)=SSvalues_AggVars_raw(SteadyStateDist, PolicyIndexes, SSvaluesFn, n_d, n_a, N_s, d_grid, a_grid,s_grid,pi_s,p); %the two zeros represent the d variables
%end
    

end