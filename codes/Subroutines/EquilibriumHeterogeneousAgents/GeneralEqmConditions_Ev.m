function  [GeneralEqmConditionsVec]=GeneralEqmConditions_Ev(SSvalues_AggVars,p, GeneralEqmEqns, Parameters, GeneralEqmEqnParamNames)

GeneralEqmConditionsVec=ones(1,length(GeneralEqmEqns))*Inf;
for i=1:length(GeneralEqmEqns)
    if isempty(GeneralEqmEqnParamNames(i).Names)  % Check for 'GeneralEqmEqnParamNames(i).Names={}'
        GeneralEqmConditionsVec(i)=GeneralEqmEqns{i}(SSvalues_AggVars, p);
    else
        GeneralEqmEqnParamsVec=CreateVectorFromParams(Parameters,GeneralEqmEqnParamNames(i).Names);
        GeneralEqmEqnParamsCell=cell(length(GeneralEqmEqnParamsVec),1);
        for jj=1:length(GeneralEqmEqnParamsVec)
            GeneralEqmEqnParamsCell(jj,1)={GeneralEqmEqnParamsVec(jj)};
        end

        GeneralEqmConditionsVec(i)=GeneralEqmEqns{i}(SSvalues_AggVars, p, GeneralEqmEqnParamsCell{:});
    end
    
end

end
