function PolicyValues=PolicyInd2Val(PolicyIndexes,n_d,n_a,n_z,d_grid)

l_d=length(n_d);

N_a=prod(n_a);
N_z=prod(n_z);

PolicyIndexes=reshape(PolicyIndexes,[l_d,N_a*N_z]);
PolicyValues=zeros(l_d,N_a*N_z);

cumsum_n_d=cumsum(n_d);

temp_d_grid=d_grid(1:cumsum_n_d(1));
PolicyValues(1,:)=temp_d_grid(PolicyIndexes(1,:));
if l_d>1
    if l_d>2
        for ii=2:l_d
            temp_d_grid=d_grid(1+cumsum_n_d(ii-1):cumsum_n_d(ii));
            PolicyValues(ii,:)=temp_d_grid(PolicyIndexes(ii,:));
        end
    end
    temp_d_grid=d_grid(cumsum_n_d(end-1)+1:end);
    PolicyValues(end,:)=temp_d_grid(PolicyIndexes(end,:));
end

PolicyValues=reshape(PolicyValues,[l_d,n_a,n_z]);


end