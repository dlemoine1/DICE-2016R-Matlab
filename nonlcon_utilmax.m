function [cneq,ceq,gneq,geq] = nonlcon_utilmax( controls, tstart, horizon, Fun, Params )

% Enforces constraints and potentially returns gradients of constraints

% Last edited: July 24, 2021 by Derek Lemoine

cneq = [];
ceq = [];
gneq = [];
geq = [];


%% Obtain trajectories

if Params.transitionsasconstraints~=1
    
    controls = controls.*Params.normalization;
    [ ~, Ynet, ~, T, M, emsind, ~, ~, abatecost] = trajectory( controls, tstart, horizon, Fun, Params );    
    
else
    
    controls = reshape(controls,horizon,[]);
    controls = controls.*Params.normalization;
    
    abaterate = controls(:,Params.col_abaterate);
    savingsrate = controls(:,Params.col_savingsrate);
    if Params.consumptionascontrol
        cons = controls(:,Params.col_C);
    end
    K = controls(:,Params.col_K);
    T = controls(:,Params.col_T);
    Tocean = controls(:,Params.col_Tocean);
    M = controls(:,Params.col_M);
    if strcmp(Params.carbonmodel,'fair') % add alpha controls    
        alpha = controls(:,Params.col_alpha);
    end    
        
    % adjust exogenous variable to start at tstart
    pop = Params.pop(tstart:tstart+horizon-1,:);
    tfp = Params.tfp(tstart:tstart+horizon-1,:);
    sigma = Params.sigma(tstart:tstart+horizon-1,:);
    psi = Params.psi(tstart:tstart+horizon-1,:);
    
    % calculate other useful variables
    Ygross = Fun.Ygross(tfp,pop,K);
    emsind = Params.timestep*sigma.*(1-abaterate).*Ygross;
    ems = emsind + Params.timestep*Fun.otherems([tstart:tstart+horizon-1]');
    Ynet = Fun.Ynet(tfp,pop,K,T);
    abatecost = Fun.abatecost(psi,abaterate,tfp,pop,K);
    inv = savingsrate.*( Ynet - abatecost );    
    
end

cumulemsind = cumsum(emsind);
cumulemsind = Params.cumulems0 + [0; cumulemsind(1:end-1,:)];
% also create a variable for constraining fossil resource use - this
% variable is not decremented by negative emissions
cumulfossil = cumsum(max(0,emsind));
cumulfossil = Params.cumulems0 + [0; cumulfossil(1:end-1,:)];


%% Establish constraints

if Params.transitionsasconstraints~=1
    
    % output constraint
    cneq = abatecost - Ynet;
    
    % cumulative fossil extraction constraint
    if Params.cumulemsindlimit<Inf % save time by not calculating if =Inf
        cneq(end+1:end+size(cumulfossil,1),1) = cumulfossil - Params.cumulemsindlimit;
    end
    
    % impulse response constraint for carbon feedback parameter
    if strcmp(Params.carbonmodel,'fair')
        alpha = controls(:,1);
        otherems = [100;Fun.otherems([tstart:tstart+horizon-1-1]')]; % T x 1 vector of cumulative emissions from deforestation (Gt C); is the leading 100 Gt C meant to adjust for all emissions prior to 2015?
        ceq = Fun.IRF1(T,M,cumulemsind + cumsum(otherems)) - Fun.IRF2(alpha); % reconcile emissions-implied and alpha-implied 100-year integrated impulse response at each time period
    end

else
    
    % set up constraints
    constr = 0;
    if Params.consumptionascontrol
        constr = constr+1;
        ceq(:,constr) = Ynet - cons - abatecost - inv; % output constraint
        if Params.scaleconstraints
            ceq(:,constr) = ceq(:,constr)./C;
        end
    end
    constr = constr+1;
    ceq(:,constr) = [ K(1,1) - Params.K0; K(2:end,1) - Fun.Knext(K(1:end-1,1),inv(1:end-1,1)) ]; % capital transition
    if Params.scaleconstraints
        ceq(:,constr) = ceq(:,constr)./K;
    end
    constr = constr+1;
    ceq(:,constr) = [ T(1,1) - Params.T0; T(2:end,1) - Fun.Tnext(T(1:end-1,1),Tocean(1:end-1,1),M(2:end,:),[tstart+1:horizon+tstart-1]') ]; % temperature transition
    if Params.scaleconstraints
        ceq(:,constr) = ceq(:,constr)./T;
    end
    constr = constr+1;
    ceq(:,constr) = [ Tocean(1,1) - Params.Tocean0; Tocean(2:end,1) - Fun.Toceannext(T(1:end-1,1),Tocean(1:end-1,1)) ]; % ocean temperature transition
    if Params.scaleconstraints
        ceq(:,constr) = ceq(:,constr)./Tocean;
    end
    constr = constr+1;
    if ~strcmp(Params.carbonmodel,'fair')
        ceq(:,constr:constr-1+size(M,2)) = [ M(1,:) - Params.M0; M(2:end,:) - Fun.Mnext(M(1:end-1,:),ems(1:end-1,1)) ];
        if Params.scaleconstraints
            ceq(:,constr:constr-1+size(M,2)) = ceq(:,constr:constr-1+size(M,2))./M;
        end
    else
        ceq(:,constr:constr-1+size(M,2)) = [ M(1,:) - Params.M0; M(2:end,:) - Fun.Mnext(M(1:end-1,:),ems(1:end-1,1),alpha(1:end-1,1)) ];
        if Params.scaleconstraints
            ceq(:,constr:constr-1+size(M,2)) = ceq(:,constr:constr-1+size(M,2))./M;
        end
        otherems = [100;Fun.otherems([tstart:tstart+horizon-1-1]')]; % T x 1 vector of cumulative emissions from deforestation (Gt C); is the leading 100 Gt C meant to adjust for all emissions prior to 2015?
        ceq(:,end+1) = Fun.IRF1(T,M,cumulemsind + cumsum(otherems)) - Fun.IRF2(alpha); % reconcile emissions-implied and alpha-implied 100-year integrated impulse response at each time period
    end
        
    % gradient of constraints: column i corresponds to row i of ceq
    geq = sparse(zeros(size(controls,1)*size(controls,2),size(ceq,1)*size(ceq,2)));
    
    % will construct controls-per-period x constraints matrices of time t and
    % time t-1 effects for each constraint, which will then use when
    % looping through time to construct the full gradient matrix
    g_thisperiod = transpose(zeros(size(controls,2),numel(ceq))); % constraints x controls-per-period; will re-transpose later
    g_lastperiod = transpose(zeros(size(controls,2),numel(ceq))); % constraints x controls-per-period; will re-transpose later            
    constr = 0;
    % gradient of output constraint
    if Params.consumptionascontrol        
        constr = constr+1;
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_K) = (1-savingsrate).*( Fun.dYnet_dK(tfp,pop,K,T) - Fun.dabatecost_dK(psi,abaterate,tfp,pop,K) );
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_T) = (1-savingsrate).*Fun.dYnet_dT(tfp,pop,K,T);
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_abaterate) = (1-savingsrate).*( - Fun.dabatecost_dabaterate(psi,abaterate,tfp,pop,K) );
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_savingsrate) = -( Ynet - abatecost );
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_C) = -1;
        if Params.scaleconstraints
            g_thisperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_thisperiod([1:horizon]+(constr-1)*horizon,:),C);
            g_lastperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_lastperiod([1:horizon]+(constr-1)*horizon,:),C);
            g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_C) = g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_C) ...
                - ceq(:,constr)./C;
        end
    end
    % gradient of capital transition
    constr = constr+1;
    g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_K) = 1; % effect on contemporaneous capital    
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_K) = -Fun.dKnext_dK(K(1:end-1,1),inv(1:end-1,1)) ...
        - Fun.dKnext_dinv(K(1:end-1,1),inv(1:end-1,1)).*savingsrate(1:end-1,1).*( Fun.dYnet_dK(tfp(1:end-1,1),pop(1:end-1,1),K(1:end-1,1),T(1:end-1,1)) - Fun.dabatecost_dK(psi(1:end-1,1),abaterate(1:end-1,1),tfp(1:end-1,1),pop(1:end-1,1),K(1:end-1,1)) ); % effect on next period's capital
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_savingsrate) = -Fun.dKnext_dinv(K(1:end-1,1),inv(1:end-1,1)).*( Ynet(1:end-1,1) - abatecost(1:end-1,1) ); % effect on next period's capital
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_T) = -Fun.dKnext_dinv(K(1:end-1,1),inv(1:end-1,1)).*savingsrate(1:end-1,1).*Fun.dYnet_dT(tfp(1:end-1,1),pop(1:end-1,1),K(1:end-1,1),T(1:end-1,1)); % effect on next period's capital
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_abaterate) = -Fun.dKnext_dinv(K(1:end-1,1),inv(1:end-1,1)).*savingsrate(1:end-1,1).*(-1).*Fun.dabatecost_dabaterate(psi(1:end-1,1),abaterate(1:end-1,1),tfp(1:end-1,1),pop(1:end-1,1),K(1:end-1,1)); % effect on next period's capital
    if Params.scaleconstraints
        g_thisperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_thisperiod([1:horizon]+(constr-1)*horizon,:),K);
        g_lastperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_lastperiod([1:horizon]+(constr-1)*horizon,:),K);
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_K) = g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_K) ...
            - ceq(:,constr)./K;
    end
    % gradient of temperature transition
    constr = constr+1;
    g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_T) = 1; % effect on contemporaneous temperature
    dTnext_dM = [zeros(1,size(M,2)); Fun.dTnext_dM(T(1:horizon-1,1),Tocean(1:horizon-1,1),M(2:horizon,:))]; % TxSinks
    for index_sink=1:size(M,2)
        g_thisperiod([2:horizon]+(constr-1)*horizon,Params.col_M(index_sink)) = -dTnext_dM(2:end,index_sink); % effect of reservoir index_sink on contemporaneous temperature; skip first period because no transition to that temperature
    end
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_T) = -Fun.dTnext_dT(T(1:horizon-1,1),Tocean(1:horizon-1,1),M(2:horizon,:)); % effect on next period's temperature
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_Tocean) = -Fun.dTnext_dTocean(T(1:horizon-1,1),Tocean(1:horizon-1,1),M(2:horizon,:)); % effect on next period's temperature    
    if Params.scaleconstraints
        g_thisperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_thisperiod([1:horizon]+(constr-1)*horizon,:),T);
        g_lastperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_lastperiod([1:horizon]+(constr-1)*horizon,:),T);
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_T) = g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_T) ...
            - ceq(:,constr)./T;
    end
    % gradient of ocean temperature transition
    constr = constr+1;
    g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_Tocean) = 1; % effect on contemporaneous temperature
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_Tocean) = -Fun.dToceannext_dTocean(T(1:horizon-1,1),Tocean(1:horizon-1,1)); % effect on next period's temperature
    g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_T) = -Fun.dToceannext_dT(T(1:horizon-1,1),Tocean(1:horizon-1,1)); % effect on next period's temperature
    if Params.scaleconstraints
        g_thisperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_thisperiod([1:horizon]+(constr-1)*horizon,:),Tocean);
        g_lastperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_lastperiod([1:horizon]+(constr-1)*horizon,:),Tocean);
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_Tocean) = g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_Tocean) ...
            - ceq(:,constr)./Tocean;
    end
    % gradient of carbon transition
    if ~strcmp(Params.carbonmodel,'fair')
        
        for index_nextsink=1:size(M,2) % looping over reservoirs are differentiating
            constr = constr + 1;
            dMnext_dM = Fun.dMnext_dM(M(1:horizon-1,:),ems(1:horizon-1,1),index_nextsink); % TxSinks
            g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) = 1; % effect on contemporaneous reservoir index_nextsink
            for index_fromsink=1:size(M,2) % looping over reservoirs are differentiating with respect to
                g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_M(index_fromsink)) = -dMnext_dM(:,index_fromsink); % effect on next period's reservoir index_nextsink
            end
            dMnext_dems = Fun.dMnext_dems(M(1:horizon-1,:),ems(1:horizon-1,1),index_nextsink);
            g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_K) = -dMnext_dems.*( Params.timestep*sigma(1:horizon-1,1).*(1-abaterate(1:horizon-1,1)).*Fun.dYgross_dK(tfp(1:horizon-1,1),pop(1:horizon-1,1),K(1:horizon-1,1))  ); % effect on next period's reservoir index_nextsink
            g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_abaterate) = -dMnext_dems.*( Params.timestep*sigma(1:horizon-1,1).*(-1).*Ygross(1:horizon-1,1)  ); % effect on next period's reservoir index_nextsink
            
            if Params.scaleconstraints
                g_thisperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_thisperiod([1:horizon]+(constr-1)*horizon,:),M(:,index_nextsink));
                g_lastperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_lastperiod([1:horizon]+(constr-1)*horizon,:),M(:,index_nextsink));
                g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) = g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) ...
                    - ceq(:,constr)./M(:,index_nextsink);
            end
            
        end        
        
    else
        
        % here dMnext_dM is Tx4, where column i gives effect on sink i and
        % we recall that sinks other than i cnnot affect sink i        
        dMnext_dM = Fun.dMnext_dM(M(1:horizon-1,:),ems(1:horizon-1,1),alpha(1:horizon-1,1)); % Tx4
        dMnext_dems = Fun.dMnext_dems(M(1:horizon-1,:),ems(1:horizon-1,1),alpha(1:horizon-1,1)); % Tx4
        dMnext_dalpha = Fun.dMnext_dalpha(M(1:horizon-1,:),ems(1:horizon-1,1),alpha(1:horizon-1,1)); % Tx4
        
        for index_nextsink=1:size(M,2) % looping over reservoirs are differentiating
            constr = constr + 1;
            g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) = 1; % effect on contemporaneous reservoir index_nextsink
            g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) = -dMnext_dM(:,index_nextsink); % effect on next period's reservoir index_nextsink
            g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_K) = -dMnext_dems(:,index_nextsink).*( Params.timestep*sigma(1:horizon-1,1).*(1-abaterate(1:horizon-1,1)).*Fun.dYgross_dK(tfp(1:horizon-1,1),pop(1:horizon-1,1),K(1:horizon-1,1))  ); % effect on next period's reservoir index_nextsink
            g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_abaterate) = -dMnext_dems(:,index_nextsink).*( Params.timestep*sigma(1:horizon-1,1).*(-1).*Ygross(1:horizon-1,1) ); % effect on next period's reservoir index_nextsink
            g_lastperiod([2:horizon]+(constr-1)*horizon,Params.col_alpha) = -dMnext_dalpha(:,index_nextsink); % effect on next period's reservoir index_nextsink
            
            if Params.scaleconstraints
                g_thisperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_thisperiod([1:horizon]+(constr-1)*horizon,:),M(:,index_nextsink));
                g_lastperiod([1:horizon]+(constr-1)*horizon,:) = bsxfun(@rdivide,g_lastperiod([1:horizon]+(constr-1)*horizon,:),M(:,index_nextsink));
                g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) = g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_nextsink)) ...
                    - ceq(:,constr)./M(:,index_nextsink);
            end
            
        end
        
        % equations implicitly defining alpha
        constr = constr + 1;        
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_T) = Fun.dIRF1_dT(T,M,cumulemsind + cumsum(otherems));
        dIRF1_dM = Fun.dIRF1_dM(T,M,cumulemsind + cumsum(otherems)); % Tx4
        for index_fromsink=1:size(M,2)
            g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_M(index_fromsink)) = dIRF1_dM(:,index_fromsink);
        end
        g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_alpha) = -Fun.dIRF2_dalpha(alpha);
        % will do these two below so can account for cumulative emissions,
        % which mean that all earlier controls enter
        % g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_K) = Fun.dIRF1_dems(T,M,cumulemsind + cumsum(otherems)).*( Params.timestep*sigma.*(1-abaterate).*Fun.dYgross_dK(tfp,pop,K)  );
        % g_thisperiod([1:horizon]+(constr-1)*horizon,Params.col_abaterate) = Fun.dIRF1_dems(T,M,cumulemsind + cumsum(otherems)).*( Params.timestep*sigma.*(-1).*Ygross  );        
        
    end
    
    ceq = ceq(:); % stack in one column vector
    g_thisperiod = g_thisperiod'; % controls-per-period x constraints
    g_lastperiod = g_lastperiod'; % controls-per-period x constraints
    % build out full gradient matrix
    for t=1:horizon
        g_thisperiod(:,t+[0:(length(ceq)/horizon-1)]*horizon) = bsxfun(@times,g_thisperiod(:,t+[0:(length(ceq)/horizon-1)]*horizon),Params.normalization(t,:)');
        geq(t+[0:(size(controls,2)-1)]*horizon, t+[0:(length(ceq)/horizon-1)]*horizon) = g_thisperiod(:,t+[0:(length(ceq)/horizon-1)]*horizon);
        if t~=1
            g_lastperiod(:,t+[0:(length(ceq)/horizon-1)]*horizon) = bsxfun(@times,g_lastperiod(:,t+[0:(length(ceq)/horizon-1)]*horizon),Params.normalization(t-1,:)');
            geq(t-1+[0:(size(controls,2)-1)]*horizon, t+[0:(length(ceq)/horizon-1)]*horizon) = g_lastperiod(:,t+[0:(length(ceq)/horizon-1)]*horizon); % note the t-1 in the row of geq: are identifying effects of t-1 controls
        end
        if strcmp(Params.carbonmodel,'fair') % account for cumulative emissions in condition defining alpha
            cum_otherems = cumsum(otherems(1:t,1));
            % requires constr not to have changed since alpha:
            geq([1:t-1]+(Params.col_K-1)*horizon, t+(constr-1)*horizon) = Fun.dIRF1_dems(T(t),M(t,:),cumulemsind(t) + cum_otherems(t)).*( Params.timestep*sigma(1:t-1,1).*(1-abaterate(1:t-1,1)).*Fun.dYgross_dK(tfp(1:t-1,1),pop(1:t-1,1),K(1:t-1,1))  ); % earlier periods' capital's effect on time t IRF1
            geq([1:t-1]+(Params.col_K-1)*horizon, t+(constr-1)*horizon) = geq([1:t-1]+(Params.col_K-1)*horizon, t+(constr-1)*horizon).*Params.normalization(1:t-1,Params.col_K);
            geq([1:t-1]+(Params.col_abaterate-1)*horizon, t+(constr-1)*horizon) = Fun.dIRF1_dems(T(t),M(t,:),cumulemsind(t) + cumsum(otherems(t))).*( Params.timestep*sigma(1:t-1,1).*(-1).*Ygross(1:t-1,1)  ); % earlier periods' abatement's effect on time t IRF1
            geq([1:t-1]+(Params.col_abaterate-1)*horizon, t+(constr-1)*horizon) = geq([1:t-1]+(Params.col_abaterate-1)*horizon, t+(constr-1)*horizon).*Params.normalization(1:t-1,Params.col_abaterate);
        end
    end    
    
    % fix certain periods' savings
    if Params.fixsavings>=0
        ceq(end+1:end+horizon,1) = savingsrate - Params.fixsavings;
        geq([1:horizon]+(Params.col_savingsrate-1)*horizon,end+1:end+horizon) = bsxfun(@times,eye(horizon),Params.normalization(1:horizon,Params.col_savingsrate));
    elseif Params.dicelrsavings==1 && Params.dicelrsavings_firstper <= size(savingsrate,1) + tstart - 1
        time_begin_fixing_savings = Params.dicelrsavings_firstper-tstart+1;
        ceq(end+1:end+horizon-time_begin_fixing_savings+1,1) = savingsrate(time_begin_fixing_savings:horizon,1) - Params.lrsavingsrate;
        geq([time_begin_fixing_savings:horizon]+(Params.col_savingsrate-1)*horizon,end+1:end+horizon-time_begin_fixing_savings+1) = bsxfun(@times,eye(horizon-time_begin_fixing_savings+1),Params.normalization([time_begin_fixing_savings:horizon],Params.col_savingsrate));
    end
    
    % cumulative fossil extraction constraint
    if Params.cumulemsindlimit<Inf % save time by not calculating if =Inf
        cneq = cumulfossil - Params.cumulemsindlimit;
        gneq = sparse(zeros(size(controls,1)*size(controls,2),size(cneq,1)));
        for t=2:horizon
            gneq([1:t-1]+(Params.col_K-1)*horizon, t) = ( Params.timestep*sigma(1:t-1,1).*(1-min(1,abaterate(1:t-1,1))).*Fun.dYgross_dK(tfp(1:t-1,1),pop(1:t-1,1),K(1:t-1,1))  ); % earlier periods' capital's effect on time t cumulative emissions
            gneq([1:t-1]+(Params.col_K-1)*horizon, t) = gneq([1:t-1]+(Params.col_K-1)*horizon, t).*Params.normalization(1:t-1,Params.col_K);
            gneq([1:t-1]+(Params.col_abaterate-1)*horizon, t) = (abaterate(1:t-1,1) < 1).*( Params.timestep*sigma(1:t-1,1).*(-1).*Ygross(1:t-1,1)  ); % earlier periods' abatement's effect on time t cumulative emissions
            gneq([1:t-1]+(Params.col_abaterate-1)*horizon, t) = gneq([1:t-1]+(Params.col_abaterate-1)*horizon, t).*Params.normalization(1:t-1,Params.col_abaterate);
        end
    end
    
    if Params.screenreport==1
        disp(['Max abs constraint deviation is ' num2str(max(abs(ceq)))]);
    end
    
end

end