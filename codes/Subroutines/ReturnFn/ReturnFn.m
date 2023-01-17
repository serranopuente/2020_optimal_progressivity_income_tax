function F=ReturnFn(lab,anext,a,z,r,gamma,varphi,chi,elle,alpha,delta,e1,e2,e3,e4,omega,lambda,tau,kappa)
% This function evaluates the utility for each state.
% Calculate the wage (from first order conditions of the firm's problem)
w=(1-alpha)*(((r+delta)/(alpha))^(alpha/(alpha-1)));

% Determine values of e (endowment of efficiency labor units) and omega_e
% (transfer to retirees) as functions of 'age' and 'labor efficiency units'
e=0;
omega_s=0;
if z==1
    e=e1;
elseif z==2
    e=e2;
elseif z==3
    e=e3;
elseif z==4
    e=e4;
else
    omega_s=omega;
end

% Calculate income
y=a*r+e*lab*w+omega_s;
% Tax on income
tau_y=(y-lambda*y^(1-tau))+kappa*y;

F=-Inf;
% Calculate consumption including income tax - augmented HSV (2017) tax function
c=a+y-tau_y-anext;

% Calculate utility or return
if c>0
    if gamma==1
        F=log(c)+chi*((elle-lab)^(1-varphi))/(1-varphi);
    else
        F=(c^(1-gamma))/(1-gamma)+chi*((elle-lab)^(1-varphi))/(1-varphi);
    end
end

end