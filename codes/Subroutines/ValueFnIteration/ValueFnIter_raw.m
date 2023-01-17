function [VKron,PolicyIndexesKron]=ValueFnIter_raw(VKron, n_d, n_a, n_z, Gamma, beta, Fmatrix, Phi_aprime,Tolerance,Howards,Howards2)

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

PolicyIndexesKron=zeros(N_a,N_z); % Indexes the optimal choice for d (labor choice and asset position for next period) given rest of dimensions a,z

currdist=Inf;
tempcounter=1; 

 while currdist>Tolerance
        
    VKronold=VKron;

    for z_c=1:N_z
        % Calculate the conditional expectation (except beta)
        RHSpart2=zeros(N_d,1);
        for zprime_c=1:N_z
            if Gamma(z_c,zprime_c)~=0 % Multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
                for d_c=1:N_d
                    RHSpart2(d_c)=RHSpart2(d_c)+VKronold(Phi_aprime(d_c,z_c,zprime_c),zprime_c)*Gamma(z_c,zprime_c);
                end
            end
        end

        for a_c=1:N_a
            entireRHS=Fmatrix(:,a_c,z_c)+beta*RHSpart2; % d by 1

            % Then maximizing d indexes
            [VKron(a_c,z_c),PolicyIndexesKron(a_c,z_c)]=max(entireRHS,[],1);
        end
    end

    VKrondist=reshape(VKron-VKronold,[N_a*N_z,1]); VKrondist(isnan(VKrondist))=0;
    currdist=max(abs(VKrondist));

    if isfinite(currdist) && tempcounter<Howards2  % Use Howards Policy Fn Iteration Improvement
        for Howards_counter=1:Howards
            VKrontemp=VKron;
            for z_c=1:N_z
                for a_c=1:N_a
                    temp=0;
                    for zprime_c=1:N_z
                        if Gamma(z_c,zprime_c)~=0 % Multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
                            temp=temp+VKrontemp(Phi_aprime(PolicyIndexesKron(a_c,z_c),z_c,zprime_c),zprime_c)*Gamma(z_c,zprime_c);
                        end
                    end
                    VKron(a_c,z_c)=Fmatrix(PolicyIndexesKron(a_c,z_c),a_c,z_c)+beta*temp;
                end
            end
        end
    end

    if rem(tempcounter,100)==0
        disp(tempcounter)
        disp(currdist)
    end

    tempcounter=tempcounter+1;
 end    
end