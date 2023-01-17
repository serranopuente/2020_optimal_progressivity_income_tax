function Gini=Gini_from_LorenzCurve(LorenzCurve)
% Takes in a LorenzCurve with N evenly spaced point points from 1/N up to 1 (in practice my
% codes default to using N=100 when creating the LorenzCurve). Returns the Gini
% coefficient.

N=length(LorenzCurve);

% Use the Gini=A/(A+B)=2*A formulation for Gini coefficent (see wikipedia).
A=0;
A=sum((1:1:N)/N-reshape(LorenzCurve,[1,N]))/N;
Gini=2*A;

end
