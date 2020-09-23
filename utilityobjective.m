function [Welfare, grad] = utilityobjective( controls, tstart, horizon, Fun, Params )

% Returns welfare (actually, negative of welfare, so can be objective of a
% minimizer).  Also can return gradient.

% Last edited: July 13, 2020 by Derek Lemoine

global iteration_number

grad = [];

iteration_number = iteration_number + 1;

totalcontrols = numel(controls);

% exogenous variables
tfp = Params.tfp([tstart:horizon+tstart-1]');
pop = Params.pop([tstart:horizon+tstart-1]');
psi = Params.psi([tstart:horizon+tstart-1]');
sigma = Params.sigma([tstart:horizon+tstart-1]');
discfactor = Params.discfactor([1:horizon]');

% Obtain consumption trajectory
if Params.transitionsasconstraints~=1
    controls = controls.*Params.normalization;
    [C, ~, K, T] = trajectory( controls, tstart, horizon, Fun, Params );
else
    controls = reshape(controls,horizon,[]);
    controls = controls.*Params.normalization;
    if Params.consumptionascontrol
        C = controls(:,Params.col_C);
    else
        abaterate = controls(:,Params.col_abaterate);
        savingsrate = controls(:,Params.col_savingsrate);
        K = controls(:,Params.col_K);
        T = controls(:,Params.col_T);
                
        Ynet = Fun.Ynet(tfp,pop,K,T);
        abatecost = Fun.abatecost(psi,abaterate,tfp,pop,K);
        inv = savingsrate.*( Ynet - abatecost );
        
        C = Ynet - abatecost - inv; % output constraint
    end
end

% utility
if min(C)>0
    Welfare = -sum(discfactor.*pop.*Fun.utility(C,pop)); % convert to present value, for first period; sum over time; negate for minimization
    % gradient of utility
    if Params.transitionsasconstraints==1
        if Params.consumptionascontrol
            grad = zeros(size(controls,1),size(controls,2));
            grad(:,Params.col_C) = -discfactor.*pop.*Fun.dutility_dC(C,pop);            
        else
            grad = zeros(size(controls,1),size(controls,2));
            
            grad(:,Params.col_K) = (1-savingsrate).*( Fun.dYnet_dK(tfp,pop,K,T) - Fun.dabatecost_dK(psi,abaterate,tfp,pop,K) );
            grad(:,Params.col_T) = (1-savingsrate).*Fun.dYnet_dT(tfp,pop,K,T);
            grad(:,Params.col_abaterate) = (1-savingsrate).*( - Fun.dabatecost_dabaterate(psi,abaterate,tfp,pop,K) );
            grad(:,Params.col_savingsrate) = -( Ynet - abatecost );
            
            % turn into (negative of) discounted marginal utility
            grad(:,:) = bsxfun(@times,-discfactor.*pop.*Fun.dutility_dC(C,pop),grad(:,:));                        
        end
        
        % adjust for normalization
        grad = grad.*Params.normalization;
            
        % stack all periods of control 1, then control 2, and so on
        grad = grad(:);
            
    end
else
    Welfare = NaN;
    if Params.transitionsasconstraints==1
        grad = NaN*ones(totalcontrols,1);
    end
end


% Save first iteration and then every 10,000th iteration
if ~isempty(iteration_number) && ( mod(iteration_number,1e4)==0 || iteration_number==1 )    
    % Create text to write
    iteration_string = cell(3,1);
    iteration_string(1,1) = {['Iteration: ' num2str(iteration_number) ', Value: ' num2str(-Welfare) '\n']};
    iteration_string(2,1) = {['Controls: ' mat2str(controls) '\n']};
    iteration_string(3,1) = {[' \n']}; % end with blank line    
    
    % Write text
    [fid,msg] = fopen([Params.savedir 'iteration_report.txt'],'at');
    if fid == -1
        error(msg);
    end
    for index_write=1:size(iteration_string,1)
        fprintf(fid, '%s \n',iteration_string{index_write,1});
    end
    fclose(fid);
end


% Reporting to screen:
if Params.screenreport==1 && Params.dohpc ~= 1
    disp(['Iteration ' num2str(iteration_number)]);
    if Params.transitionsasconstraints~=1
        disp(['Abatement rate controls are ' mat2str(controls(1:horizon,1))]);
    else
        disp(['Abatement rate controls are ' mat2str(controls(:,Params.col_abaterate))]);
    end
    disp(['Value is ' num2str(-Welfare), ' with minimum consumption of ' num2str(min(C))]);
end



end

