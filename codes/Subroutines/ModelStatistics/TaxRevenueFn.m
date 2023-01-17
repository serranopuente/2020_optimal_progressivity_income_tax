function [IncomeTaxRevenue]=TaxRevenueFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa)

% lab is 'h' (hours worked)
% a is asset position or capital holding in current period
% z determines age & retirement

% Calculate Wage (from FOC of the firm's maximization problem)

w=(1-alpha)*(((r+delta)/(alpha))^(alpha/(alpha-1)));

% Calculate Income
    % If retired
    e=0;
    y=a*r+e*lab*w+omega;
    %If working
    if z<=J 
    e=e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4);
    y=a*r+e*lab*w;
    end

% Calculate Income Tax Revenue - augmented HSV (2017) tax function
IncomeTaxRevenue=(y-lambda*y^(1-tau))+kappa*y;

end