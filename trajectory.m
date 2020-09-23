function [ C, Ynet, K, T, M, emsind, Tocean, Ygross, abatecost] = trajectory( controls, tstart, horizon, Fun, Params )
 
% Simulates from tstart for horizon periods, given vector of controls.  Returns trajectories of interest.

% Last edited: July 6, 2020 by Derek Lemoine

% Capture inputs
if strcmp(Params.carbonmodel,'fair')
    alpha = controls(1:horizon,1);
    controls(1:horizon,:) = [];
end
abaterate = controls(1:horizon,1);
if Params.fixsavings<0
    savingsrate = controls(horizon+1:end,1);
else
    savingsrate = Params.fixsavings*ones(size(abaterate,1),1);
end

% Initialize variables
C = zeros(horizon,1);
Ynet = zeros(horizon,1);
K = zeros(horizon,1);
T = zeros(horizon,1);
M = zeros(horizon,length(Params.M0));
emsind = zeros(horizon,1);
Tocean = zeros(horizon,1);
Ygross = zeros(horizon,1);
abatecost = zeros(horizon,1);
inv = zeros(horizon,1);
ems = zeros(horizon,1);

% Loop through time
for t=1:horizon
                
    periods_from_start = t + tstart - 1; % current period, in absolute terms
    
    % Update time t variables
    if t==1 % get starting variables
        
        K(t,1) = Params.K0;
        M(t,:) = Params.M0;
        T(t,1) = Params.T0;
        Tocean(t,1) = Params.Tocean0;
        
    else % transitions
        
        K(t,1) = Fun.Knext(K(t-1,1),inv(t-1,1));
        if ~strcmp(Params.carbonmodel,'fair')
            M(t,:) = Fun.Mnext(M(t-1,:),ems(t-1,1));
        else
            M(t,:) = Fun.Mnext(M(t-1,:),ems(t-1,1),alpha(t-1,1));
        end
        T(t,1) = Fun.Tnext(T(t-1,1),Tocean(t-1,1),M(t,:),periods_from_start);
        Tocean(t,1) = Fun.Toceannext(T(t-1,1),Tocean(t-1,1));
        
    end
    
    % Determine emissions, output, and consumption
    Ygross(t,1) = Fun.Ygross(Params.tfp(periods_from_start,1),Params.pop(periods_from_start,1),K(t,1));
    emsind(t,1) = Params.timestep*Params.sigma(periods_from_start,1)*(1-abaterate(t,1))*Ygross(t,1);
    ems(t,1) = emsind(t,1) + Params.timestep*Fun.otherems(periods_from_start);    
    Ynet(t,1) = Fun.Ynet(Params.tfp(periods_from_start,1),Params.pop(periods_from_start,1),K(t,1),T(t,1));
    abatecost(t,1) = Fun.abatecost(Params.psi(periods_from_start,1),abaterate(t,1),Params.tfp(periods_from_start,1),Params.pop(periods_from_start,1),K(t,1));
    if Params.dicelrsavings==1 && periods_from_start >= Params.dicelrsavings_firstper % fixed long-run savings rate for final periods; essentially transversality condition
        inv(t,1) = Params.lrsavingsrate*( Ynet(t,1) - abatecost(t,1) );
    else
        inv(t,1) = savingsrate(t,1)*( Ynet(t,1) - abatecost(t,1) );
    end
    C(t,1) = Ynet(t,1) - inv(t,1) - Fun.abatecost(Params.psi(periods_from_start,1),abaterate(t,1),Params.tfp(periods_from_start,1),Params.pop(periods_from_start,1),K(t,1));
            
end % end looping through time


end

