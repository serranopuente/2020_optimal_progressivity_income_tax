function Fmatrix=CreateReturnFnMatrix(ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid, ReturnFnParamsVec)

disp(' ')
disp("Creating Return Matrix...")

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

d_gridvals=zeros(N_d,length(n_d));
for i1=1:N_d
    sub=zeros(1,length(n_d));
    sub(1)=rem(i1-1,n_d(1))+1;
    for ii=2:length(n_d)-1
        sub(ii)=rem(ceil(i1/prod(n_d(1:ii-1)))-1,n_d(ii))+1;
    end
    sub(length(n_d))=ceil(i1/prod(n_d(1:length(n_d)-1)));
    
    if length(n_d)>1
        sub=sub+[0,cumsum(n_d(1:end-1))];
    end
    d_gridvals(i1,:)=d_grid(sub);
end

z_gridvals=zeros(N_z,length(n_z));
for i2=1:N_z
    sub=zeros(1,length(n_z));
    sub(1)=rem(i2-1,n_z(1))+1;
    for ii=2:length(n_z)-1
        sub(ii)=rem(ceil(i2/prod(n_z(1:ii-1)))-1,n_z(ii))+1;
    end
    sub(length(n_z))=ceil(i2/prod(n_z(1:length(n_z)-1)));
    
    if length(n_z)>1
        sub=sub+[0,cumsum(n_z(1:end-1))];
    end
    z_gridvals(i2,:)=z_grid(sub);
end

a_gridvals=zeros(N_a,length(n_a));
for i3=1:N_a
    sub=zeros(1,length(n_a));
    sub(1)=rem(i3-1,n_a(1))+1;
    for ii=2:length(n_a)-1
        sub(ii)=rem(ceil(i3/prod(n_a(1:ii-1)))-1,n_a(ii))+1;
    end
    sub(length(n_a))=ceil(i3/prod(n_a(1:length(n_a)-1)));
    
    if length(n_a)>1
        sub=sub+[0,cumsum(n_a(1:end-1))];
    end
    a_gridvals(i3,:)=a_grid(sub);
end

ReturnFnParamsCell=num2cell(ReturnFnParamsVec);

Fmatrix=zeros(N_d,N_a,N_z);
    parfor i3=1:N_z
        z_gridvals_temp=z_gridvals(i3,:);
        Fmatrix_z=zeros(N_d,N_a);
        for i1=1:N_d
            for i2=1:N_a
                Fmatrix_z(i1,i2)=ReturnFn(d_gridvals(i1,1),d_gridvals(i1,2),a_gridvals(i2,1),z_gridvals_temp,ReturnFnParamsCell{:});
            end
        end
        Fmatrix(:,:,i3)=Fmatrix_z;
    end
    
disp(' ')

end


