function c=ConsumptionFn(lab,anext,a,z,J,r,alpha,delta,omega,e1,e2,e3,e4,lambda,tau,kappa)

w=(1-alpha)*(((r+delta)/(alpha))^(alpha/(alpha-1)));

% Income
% If retired
e=0;
y=a*r+e*lab*w+omega;

if z<=J %If working
    e=e1*(z==1)+e2*(z==2)+e3*(z==3)+e4*(z==4);
    y=a*r+e*lab*w;
end

c=a+lambda*y^(1-tau)+(1-kappa)*y-anext;

end