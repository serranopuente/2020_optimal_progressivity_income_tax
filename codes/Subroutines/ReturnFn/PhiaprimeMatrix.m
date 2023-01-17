function Phi_aprimeKron=PhiaprimeMatrix(n_d,n_z)
% Next period's assets are a function of endogenous states (current labor supply
% and current period's assets), the current period's age today and the next period's age)

N_d=prod(n_d);  % Get the number of possible combinations of endogenous variables

Phi_aprimeKron=zeros([N_d,n_z,n_z]); % Generate matrix to store the results

% Calculate next period's assets for each possible combination of values of the endogenous variables
for i=1:N_d 
    % Next period's assets is what the chooses
    d_sub=ind2sub_homemade(n_d,i);
    d2_c=d_sub(2);
    Phi_aprimeKron_d=d2_c*ones([1,n_z,n_z]);
    % Unless he/she has to pay the inhertiance tax, in which case... (tauE
    % is state tax and zlowebar is the extemption limit)
%     if a_grid(d2_c)>=zlowerbar %if the inheritance tax is relevant
%         tau_e=tauE*(a_grid(d2_c)-zlowerbar);
%         aprime=a_grid(d2_c)-tau_e;
%         aprime_c=dsearchn(a_grid,aprime);
%         Phi_aprimeKron_d(1,J+1:2*J,1:J)=aprime_c*ones(size(Phi_aprimeKron_d(1,J+1:2*J,1:J)));
%     end
    Phi_aprimeKron(i,:,:)=Phi_aprimeKron_d;
end
end