function Dist=MethodSimMomentsFn(CalibParamsVec, CalibStats_Model_Fn, CalibStats_Data, CalibDistType, CalibWeights, NumOfCalibStats)

% Takes in current parameter values, model statistics function, data 
% statistics, a vector telling it how to calculate each of the distances 
% (ie. absolute distance, percentage dist, etc.), and a vector to give
% each of them weights.

% Outputs a vector of weighted distances (often you will just want the sum
% of this output).

% CalibParamsVec: a column vector containing the parameters to be calibrated
% CalibStats_Model_Fn: a function which takes as inputs the CalibParams, and
% gives as outputs the ModelMoments
% DataMoments: a vector containing the target data moments
% CalibStats_Data: a vector (length equal to number of model/data moments)
% containing numbers in range 1-4 that say how to calculate the distances
% (ie. aboslute distance, percentage dist., etc.)
% CalibWeights: a vector containing the weights to apply to each of the
%individual moment distances

%% Calculate CalibStats_Model from current values of the CalibParams, using ModelMomentsFn
CalibStats_Model=CalibStats_Model_Fn(CalibParamsVec);

%% Calculate the weighted distance of the calibration statistics from their targets
Dist=zeros(NumOfCalibStats,1);

for i=1:NumOfCalibStats
    if CalibDistType(i)==1
        if CalibStats_Data(i)~=0
            Dist(i)=((CalibStats_Model(i)-CalibStats_Data(i))^2)/abs(CalibStats_Data(i));
        else
            disp(sprintf('WARNING: Calibration Statistic number %d is zero valued, so have overriden corresponding choice for CalibDistType', i))
            Dist(i)=((CalibStats_Model(i)-CalibStats_Data(i))^2);
        end
    elseif CalibDistType(i)==2
        if CalibStats_Data(i)~=0
            Dist(i)=abs(CalibStats_Model(i)-CalibStats_Data(i))/abs(CalibStats_Data(i));
        else
            disp(sprintf('WARNING: Calibration Statistic number %d is zero valued, so have overriden corresponding choice for CalibDistType', i))
            Dist(i)=abs(CalibStats_Model(i)-CalibStats_Data(i));
        end
    elseif CalibDistType(i)==3
        Dist(i)=(CalibStats_Model(i)-CalibStats_Data(i))^2;
    elseif CalibDistType(i)==4
        Dist(i)=abs(CalibStats_Model(i)-CalibStats_Data(i));    
    end
end

if iscolumn(CalibWeights)
    Dist=Dist.*CalibWeights;
else
    disp(' ')
    disp('CalibWeights, in the Method of Simulated Moments function should be a column vector')
    Dist=Dist.*CalibWeights';
end
    disp(' ')
    disp('Current vector of weighted distances:')
    disp(Dist')
    fprintf('Current sum of weighted distances is: %8.4f \n', sum(Dist))

% Save the current status of the Method of Simulated Moments, this allows you to check how the estimation is progressing while it is running 
save ./Output/Calib/MethodSimMoments_Status.mat CalibParamsVec CalibStats_Data CalibStats_Model Dist

% Save sequence of convergence of algorithm
update = zeros(1,4*NumOfCalibStats+1);
update(1:NumOfCalibStats) = CalibParamsVec';
update((NumOfCalibStats+1):(2*NumOfCalibStats)) = CalibStats_Data';
update((2*NumOfCalibStats+1):(3*NumOfCalibStats)) = CalibStats_Model';
update((3*NumOfCalibStats+1):(4*NumOfCalibStats)) = Dist';
load ./Output/Calib/CalibCounter.mat CalibCounter
update(4*NumOfCalibStats+1) = CalibCounter;
filename = 'Output\Tables\Calib_algorithm.xlsx';
my_cell = sprintf( 'A%s',num2str(CalibCounter));
writematrix(update,filename,'Sheet','Sequence','Range',my_cell)

end