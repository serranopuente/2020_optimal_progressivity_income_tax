function SSvalues_MeanMedianStdDev=SSvalues_MeanMedianStdDev (StationaryDist, PolicyIndexes, FnsToEvaluateFn, Parameters,FnsToEvaluateFnParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid)
% Evaluates the mean value (weighted sum/integral), median value, and standard deviation for each element of SSvaluesFn

l_d=length(n_d);
l_a=length(n_a);
l_z=length(n_z);

N_a=prod(n_a);
N_z=prod(n_z);

SSvalues_MeanMedianStdDev=zeros(length(FnsToEvaluateFn),3); % 3 Columns: Mean, Median, and Standard Deviation

StationaryDistVec=reshape(StationaryDist,[N_a*N_z,1]);

StationaryDistVec=gather(StationaryDistVec);

SSvalues_AggVars=zeros(length(FnsToEvaluateFn),1);
d_val=zeros(l_d,1);

z_gridvals=cell(N_z,l_z);
for i1=1:N_z
    sub=zeros(1,l_z);
    sub(1)=rem(i1-1,n_z(1))+1;
    for ii=2:l_z-1
        sub(ii)=rem(ceil(i1/prod(n_z(1:ii-1)))-1,n_z(ii))+1;
    end
    sub(l_z)=ceil(i1/prod(n_z(1:l_z-1)));

    if l_z>1
        sub=sub+[0,cumsum(n_z(1:end-1))];
    end
    z_gridvals(i1,:)=num2cell(z_grid(sub));
end
a_gridvals=cell(N_a,l_a);
for i2=1:N_a
    sub=zeros(1,l_a);
    sub(1)=rem(i2-1,n_a(1))+1;
    for ii=2:l_a-1
        sub(ii)=rem(ceil(i2/prod(n_a(1:ii-1)))-1,n_a(ii))+1;
    end
    sub(l_a)=ceil(i2/prod(n_a(1:l_a-1)));

    if l_a>1
        sub=sub+[0,cumsum(n_a(1:end-1))];
    end
    a_gridvals(i2,:)=num2cell(a_grid(sub));
end

PolicyIndexes=reshape(PolicyIndexes,[size(PolicyIndexes,1),N_a,N_z]);
d_gridvals=cell(N_a*N_z,l_d);
for ii=1:N_a*N_z
    j1=rem(ii-1,N_a)+1;
    j2=ceil(ii/N_a);
    d_ind=PolicyIndexes(:,j1,j2);
    for kk1=1:l_d
        if kk1==1
            d_val(kk1)=d_grid(d_ind(kk1));
        else
            d_val(kk1)=d_grid(d_ind(kk1)+sum(n_d(1:kk1-1)));
        end
    end
    d_gridvals(ii,:)=num2cell(d_val);
end

for i=1:length(FnsToEvaluateFn)
    % Includes check for cases in which no parameters are actually required
    if isempty(FnsToEvaluateFnParamNames(i).Names) % check for 'SSvalueParamNames={}'
        Values=zeros(N_a*N_z,1);
        for ii=1:N_a*N_z
            j1=rem(ii-1,N_a)+1;
            j2=ceil(ii/N_a);
            Values(ii)=FnsToEvaluateFn{i}(d_gridvals{j1+(j2-1)*N_a,:},a_gridvals{j1,:},z_gridvals{j2,:});
        end
        % When evaluating value function (which may sometimes give -Inf
        % values) on StationaryDistVec (which at those points will be
        % 0) we get 'NaN'. Use temp as intermediate variable just eliminate those.
        temp=Values.*StationaryDistVec;
        SSvalues_AggVars(i)=sum(temp(~isnan(temp)));
    else
        SSvalueParamsCell=num2cell(CreateVectorFromParams(Parameters,FnsToEvaluateFnParamNames(i).Names));
        Values=zeros(N_a*N_z,1);
        for ii=1:N_a*N_z
            j1=rem(ii-1,N_a)+1;
            j2=ceil(ii/N_a);
            Values(ii)=FnsToEvaluateFn{i}(d_gridvals{j1+(j2-1)*N_a,:},a_gridvals{j1,:},z_gridvals{j2,:},SSvalueParamsCell{:});
        end
        % When evaluating value function (which may sometimes give -Inf
        % values) on StationaryDistVec (which at those points will be
        % 0) we get 'NaN'. Use temp as intermediate variable just eliminate those.
    end
    % Mean
    SSvalues_MeanMedianStdDev(i,1)=sum(Values.*StationaryDistVec);
    % Median
    [SortedValues,SortedValues_index] = sort(Values);
    SortedStationaryDistVec=StationaryDistVec(SortedValues_index);
    median_index=find(cumsum(SortedStationaryDistVec)>=0.5,1,'first');
    SortedValues(median_index);
    SSvalues_MeanMedianStdDev(i,2)=SortedValues(median_index);
    % Standard Deviation
    SSvalues_MeanMedianStdDev(i,3)=sqrt(sum(StationaryDistVec.*((Values-SSvalues_MeanMedianStdDev(i,1).*ones(N_a*N_z,1)).^2)));
    
end

end

