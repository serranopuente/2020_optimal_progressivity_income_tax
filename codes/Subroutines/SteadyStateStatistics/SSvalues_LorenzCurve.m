function SSvalues_LorenzCurve=SSvalues_LorenzCurve(StationaryDist, PolicyIndexes, SSvaluesFn, Parameters,SSvalueParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, npoints)
% Returns a Lorenz Curve 100-by-1 that contains all of the quantiles from 1 to 100

% Note that to unnormalize the Lorenz Curve you can just multiply it be the
% SSvalues_AggVars for the same variable. This will give you the inverse
% cdf.

if nargin<13
    npoints=100;
end

l_d=length(n_d);
l_a=length(n_a);
l_z=length(n_z);

N_a=prod(n_a);
N_z=prod(n_z);

StationaryDistVec=reshape(StationaryDist,[N_a*N_z,1]);

SSvalues_AggVars=zeros(length(SSvaluesFn),1);
SSvalues_LorenzCurve=zeros(length(SSvaluesFn),npoints);
d_val=zeros(l_d,1);
      
for i=1:length(SSvaluesFn)
    SSvalueParamsVec=CreateVectorFromParams(Parameters,SSvalueParamNames(i).Names);
    SSvalueParamsCell=num2cell(SSvalueParamsVec);
    Values=zeros(N_a,N_z);
    for j1=1:N_a
        a_ind=ind2sub_homemade([n_a],j1);
        a_val=a_grid(a_ind);
        for j2=1:N_z
            s_ind=ind2sub_homemade([n_z],j2);
            z_val=z_grid(s_ind);
            d_ind=PolicyIndexes(1:l_d,j1,j2);
            for kk1=1:l_d
                if kk1==1
                    d_val(kk1)=d_grid(d_ind(kk1));
                else
                    d_val(kk1)=d_grid(d_ind(kk1)+sum(n_d(1:kk1-1)));
                end
            end
            Values(j1,j2)=SSvaluesFn{i}(d_val(1),d_val(2),a_val,z_val,SSvalueParamsCell{:});
        end
    end
    Values=reshape(Values,[N_a*N_z,1]);

    WeightedValues=Values.*StationaryDistVec;
    SSvalues_AggVars(i)=sum(WeightedValues);


    [~,SortedValues_index] = sort(Values);

    SortedStationaryDistVec=StationaryDistVec(SortedValues_index);
    SortedWeightedValues=WeightedValues(SortedValues_index);

    CumSumSortedStationaryDistVec=cumsum(SortedStationaryDistVec);
    SSvalues_LorenzCurve(i,:)=LorenzCurve_subfunction_PreSorted(SortedWeightedValues,CumSumSortedStationaryDistVec,npoints)';
end


end

