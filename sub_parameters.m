% Parameterization for DICE-2016R

% Last edited: February 24, 2021 by Derek Lemoine


%% Parameters

Params.horizon = Params.numyears/Params.timestep; % number of periods over which will optimize policy; default: 500/Params.timestep

Params.discrate = 0.015; % annual utility discount rate; default: 0.015
Params.rra = 1.45; % coefficient of relative risk aversion; default: 1.45

Params.capshare = 0.3; %capital share in production; default: 0.3
Params.deprec = 0.1; % annual capital depreciation rate; default: 0.1
Params.lrsavingsrate = Params.capshare*(Params.deprec+0.004)/(Params.deprec+0.004*Params.rra + Params.discrate); % savings rate that will hold for final 10 periods if Params.dicelrsavings==1; avoids savings crashing as end of world nears
Params.dicelrsavings_firstper = Params.horizon-round(5*10/Params.timestep)+1; % period in which will begin fixing long-run savings rate
Params.K0 = 223; % year 2015 capital (trillion 2010$); default: 223

switch Params.damagemodel 
    case 'dice'
        Params.damcoeff = 0.00236; % damage coefficient; default: 0.00236
    case 'expert' % from Lemoine (2021) fit to Pindyck (2019) variance
        Params.damcoeff = 0.0228; % damage coefficient; default: 0.0228
    otherwise
        error('Unrecognized Params.damagemodel: %s',Params.damagemodel);
end

Params.abate_exp = 2.6; % abatement cost function exponent; default: 2.6
Params.backstop = 2016.7; % backstop cost in 2015 (2010$ per ton C); default: 2016.7
Params.gpsi = 0.025; % initial decline in backstop cost per 5 years; default: 0.025
Params.bound_abate = ones(Params.horizon,1); % upper bound on abatement rate
if Params.dicenegems==1
    Params.bound_abate(ceil(30*5/Params.timestep):Params.horizon,1) = 1.2; % later can have negative emissions
end
if Params.dicefirstperabate==1
    Params.fixabate1=0.03; % DICE-2016R fixes first-period abatement at 3%    
end

Params.L0 = 7403; % year 2015 population (millions)
Params.Linf = 11500; % asymptotic population (millions)
Params.gL = 0.134*(Params.timestep/5); % rate of approach to asymptotic population level; default: 0.134

Params.tfp0 = 5.115; % initial tfp; default: 5.115
Params.gA0 = 0.076; % initial growth rate of tfp per 5 years; default: 0.076
Params.deltaA = 0.005; % annual decline in growth rate of tfp; default: 0.005

Params.sigma0 = 0.0955; % initial emission intensity of output (Gt C per trillion 2010$); default: 0.0955
Params.gsigma0 = -0.0152; % initial annual growth rate of emission intensity; default: -0.0152
Params.deltasigma = -0.001; % annual change in growth rate of emission intensity; default: -0.001

Params.otherems0 = 0.71; % annual emissions from deforestation, in Gt C; default: 0.71
Params.gotherems = 0.115; % decline rate of deforestation emissions per 5 years; default: 0.115

Params.otherforcing0 = 0.5; % year 2015 non-CO2 forcing (W/m^2); default: 0.5
Params.otherforcing100 = 1; % year 2100 non-CO2 forcing (W/m^2); default: 1

switch Params.carbonmodel
    case 'dice' % DICE-2016R
        Params.Mmatrix = [1-0.12 0.12*588/360 0; 0.12 1-0.12*588/360-0.007 0.007*360/1720; 0 0.007 1-0.007*360/1720]; % transition matrix for carbon
        Params.Mmatrix = Params.Mmatrix^(1/5); % annualize
        Params.emsinks = [1 0 0]; % 1x3 allocation of emissions to CO2 reservoirs; default: [1 0 0]
        Params.M0 = [ 851 460 1740 ]; % 1x3; year 2015 CO2 stocks (Gt C); default: [ 851 460 1740 ]
    case 'gol' % Golosov et al. (2014), via Lemoine (2020) and with starting stocks from Dietz et al. (2020)
        Params.Mmatrix = bsxfun(@times,eye(2,2),[1 1-0.0023]); % 2x2 CO2 transition matrix; default: bsxfun(@times,eye(2,2),[1 1-0.0023])
        Params.emsinks = [0.2 0.393*(1-0.2)]; % 1x2 allocation of emissions to CO2 reservoirs; default: [0.2 0.393*(1-0.2)] - is ok that doesn't add up to 1
        Params.M0 = [ 712 159 ]; % 1x2; year 2015 CO2 stocks (Gt C); default: [ 712 159 ]
    otherwise % Joos et al. (2013), via Dietz et al. (2020); FAIR will use this too
        Params.Mmatrix = bsxfun(@times,eye(4,4),[1 0.9975 0.9730 0.7927]); % 4x4 CO2 transition matrix; default: bsxfun(@times,eye(4,4),[1 0.9975 0.9730 0.7927])
        Params.emsinks = [0.2173 0.2240 0.2824 0.2763]; % 1x4 allocation of emissions to CO2 reservoirs; default: [0.2173 0.2240 0.2824 0.2763]
        Params.M0 = [ 588+139.1 90.2 29.2 4.2 ]; % 1x4; year 2015 CO2 stocks (Gt C); default: [ 588+139.1 90.2 29.2 4.2 ]
end

switch Params.climatemodel
    case 'dice' % DICE-2016R
        Params.f2x = 3.6813; % forcing from doubled CO2 (W/m^2); default: 3.6813
        Params.lambda = Params.f2x/3.1; % forcing per degree warming (W/m^2/K); default: Params.f2x/3.1
        Params.phi1 = 0.1005;  % warming delay parameter; default: 0.1005
        Params.phi3 = 0.088; % xfer of heat from ocean to surface; default: 0.088
        Params.phi4 = 0.025; % xfer of heat from surface to ocean; default: 0.025
    otherwise % Geoffroy et al. (2013) via Lemoine (2020)
        Params.f2x = 3.503; % forcing from doubled CO2 (W/m^2); default: 3.503
        Params.lambda = 1.13; % forcing per degree warming (W/m^2/K); default: 1.13
        Params.phi1 = 0.386;  % warming delay parameter; default: 0.386
        Params.phi3 = 0.73; % xfer of heat from ocean to surface; default: 0.73
        Params.phi4 = 0.034; % xfer of heat from surface to ocean; default: 0.034
end
Params.T0 = 0.85; % year 2015 surface temperature (deg C, rel to 1900); default: 0.85
Params.Tocean0 = 0.0068; % year 2015 lower ocean temperature (deg C, rel to 1900); default: 0.0068

Params.cumulems0 = 400; % Gt C
if Params.dicecumulems==1
    Params.cumulemsindlimit = 6000; % maximum cumulative extraction of fossil fuels (Gt C)
else
    Params.cumulemsindlimit = Inf; % no limit 
end

if strcmp(Params.carbonmodel,'fair')
   Params.decay = [ 0 0.00254 0.0274 0.232342 ]; % 1x4 vector of decay in permanent, slow, medium, then fast boxes; approximately equal to -log(diag(Params.Mmatrix))
end
% terms related to carbon cycle feedbacks, from Dietz et al. (2020)
Params.R0 = 34.4; % pre-industrial iIRF (years); default: 34.4
Params.RC = 0.019; % positive feedback from cumulative carbon uptake (years/GtC); default: 0.019
Params.RT = 4.165; % positive feedback from warming (years/degC); default: 4.165


%% Functions

% Functions of time t choices
Fun.Ygross = @(tfp,L,K) tfp.*((L/1000).^(1-Params.capshare)).*(K.^Params.capshare);
Fun.Ynet = @(tfp,L,K,T) Fun.Ygross(tfp,L,K).*(1-Params.damcoeff*T.^2);
Fun.utility = @(C,L) ((C./L).^(1-Params.rra))/(1-Params.rra);
Fun.abatecost = @(psi,abaterate,tfp,L,K) psi.*(abaterate.^Params.abate_exp).*Fun.Ygross(tfp,L,K);

% Derivatives, for gradient calculations
Fun.dYgross_dK = @(tfp,L,K) Params.capshare*Fun.Ygross(tfp,L,K)./K;
Fun.dYnet_dK = @(tfp,L,K,T) Fun.Ynet(tfp,L,K,T).*Fun.dYgross_dK(tfp,L,K)./Fun.Ygross(tfp,L,K);
Fun.dYnet_dT = @(tfp,L,K,T) -2*Params.damcoeff*T.*Fun.Ygross(tfp,L,K);
Fun.dutility_dC = @(C,L) ((C./L).^(-Params.rra))./L; % derivative wrt consumption
Fun.dabatecost_dabaterate = @(psi,abaterate,tfp,L,K) Params.abate_exp*psi.*(abaterate.^(Params.abate_exp-1)).*Fun.Ygross(tfp,L,K);
Fun.dabatecost_dK = @(psi,abaterate,tfp,L,K) Fun.abatecost(psi,abaterate,tfp,L,K).*Fun.dYgross_dK(tfp,L,K)./Fun.Ygross(tfp,L,K);

% Exogenous transitions: t is number of periods, with t=1 being 2015
Fun.popnext = @(L) L.*(Params.Linf./L).^Params.gL;
Fun.gA = @(t) Params.gA0*exp(-Params.timestep*Params.deltaA*(t-1)); % deltaA is annual decline in growth rate, which matches DICE-2016R code although not parameter definition
Fun.tfpnext = @(t,tfp) tfp./((1-Fun.gA(t)).^(Params.timestep/5)); % growth rate gA is defined per timestep
Fun.gsigma = @(t) Params.gsigma0*((1+Params.deltasigma).^(Params.timestep*(t-1)));
Fun.sigmanext = @(t,sigma) sigma.*exp(Params.timestep*Fun.gsigma(t));
Fun.psi = @(t,sigma) Params.backstop*((1-Params.gpsi).^((t-1)*Params.timestep/5)).*sigma/(1000*Params.abate_exp);
Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5);
Fun.otherforcing = @(t) Params.otherforcing0 + (Params.otherforcing100 - Params.otherforcing0)*min((t-1)*Params.timestep/(17*5),1);

% Endogenous transitions
Fun.Knext = @(Know,invnow) Know*((1-Params.deprec)^Params.timestep) + Params.timestep*invnow; % capital
switch Params.carbonmodel
    case 'fair'
        % Mnow is Tx4, emsnow is Tx1, alpha is Tx1
        Fun.Mpermnext = @(Mnow,emsnow,alpha) Params.emsinks(1)*emsnow + Mnow(:,1); % box 1 - no decay; Tx1
        Fun.Mdecayingboxesnext = @(Mnow,emsnow,alpha,index) (Params.emsinks(index)./(Params.decay(index)./alpha)).*(1-exp(-Params.timestep*(Params.decay(index)./alpha))).*emsnow ...
            + exp(-Params.timestep*(Params.decay(index)./alpha)).*Mnow(:,index); % index indicates box: 2 (slow decay), 3 (medium decay), or 4 (fast decay); Tx1
        Fun.Mnext = @(Mnow,emsnow,alpha) [ Fun.Mpermnext(Mnow,emsnow,alpha) Fun.Mdecayingboxesnext(Mnow,emsnow,alpha,2) Fun.Mdecayingboxesnext(Mnow,emsnow,alpha,3) Fun.Mdecayingboxesnext(Mnow,emsnow,alpha,4)]; % CO2 reservoirs (Gt C); Tx4
    otherwise
        % Mnow is TxSinks, emsnow is Tx1, output is TxSinks
        Fun.Mnext = @(Mnow,emsnow) transpose((Params.Mmatrix^Params.timestep)*Mnow') + emsnow*Params.emsinks; % CO2 reservoirs (Gt C)
end
switch Params.carbonmodel
    case 'dice' % atmospheric CO2 is the first reservoir
        Fun.forcing = @(M,t) Params.f2x*log(M(:,1)/588)/log(2) + Fun.otherforcing(t); % forcing (W/m^2)
    otherwise % atmospheric CO2 is sum of all reservoirs
        Fun.forcing = @(M,t) Params.f2x*log(sum(M,2)/588)/log(2) + Fun.otherforcing(t); % forcing (W/m^2)
end
Fun.Tnext = @(Tnow,Toceannow,Mnext,tnext) Tnow + (Params.timestep/5)*Params.phi1*( Fun.forcing(Mnext,tnext) - Params.lambda*Tnow - Params.phi3*(Tnow - Toceannow) ); % surface temperature (deg C, wrt 1900); really should recalibrate if Params.timestep~=5
Fun.Toceannext = @(Tnow,Toceannow) Toceannow + (Params.timestep/5)*Params.phi4*( Tnow - Toceannow ); % lower ocean temperature (deg C, wrt 1900); really should recalibrate if Params.timestep~=5

% Derivatives, for gradient calculations
Fun.dKnext_dK = @(Know,invnow) ((1-Params.deprec)^Params.timestep)*ones(size(Know,1),1); % Tx1
Fun.dKnext_dinv = @(Know,invnow) Params.timestep*ones(size(invnow,1),1); % Tx1
switch Params.carbonmodel
    case 'fair'
        % Mnow is Tx4, emsnow is Tx1, alpha is Tx1
        
        % Tx1:
        Fun.dMpermnext_dM = @(Mnow,emsnow,alpha) ones(size(Mnow,1),1);
        Fun.dMpermnext_dems = @(Mnow,emsnow,alpha) Params.emsinks(1)*ones(size(emsnow,1),1);
        Fun.dMpermnext_dalpha = @(Mnow,emsnow,alpha) zeros(size(Mnow,1),1);
        
        % Tx1:
        Fun.dMdecayingboxesnext_dM = @(Mnow,emsnow,alpha,index) exp(-Params.timestep*(Params.decay(index)./alpha)).*ones(size(Mnow,1),1); % index indicates box: 2 (slow decay), 3 (medium decay), or 4 (fast decay); Tx1
        Fun.dMdecayingboxesnext_dems = @(Mnow,emsnow,alpha,index) (Params.emsinks(index)./(Params.decay(index)./alpha)).*(1-exp(-Params.timestep*(Params.decay(index)./alpha))).*ones(size(emsnow,1),1); % index indicates box: 2 (slow decay), 3 (medium decay), or 4 (fast decay); Tx1
        Fun.dMdecayingboxesnext_dalpha = @(Mnow,emsnow,alpha,index) (Params.emsinks(index)./(Params.decay(index))).*(1-exp(-Params.timestep*(Params.decay(index)./alpha))).*emsnow ...
            + (Params.emsinks(index)./(Params.decay(index)./alpha)).*(-Params.timestep*(Params.decay(index)./(alpha.^2))).*(exp(-Params.timestep*(Params.decay(index)./alpha))).*emsnow ...
            + (Params.timestep*(Params.decay(index)./(alpha.^2))).*exp(-Params.timestep*(Params.decay(index)./alpha)).*Mnow(:,index); % index indicates box: 2 (slow decay), 3 (medium decay), or 4 (fast decay); Tx1
        
        % Tx4: (use with knowledge that one stock cannot affect others)
        Fun.dMnext_dM = @(Mnow,emsnow,alpha) [ Fun.dMpermnext_dM(Mnow,emsnow,alpha) Fun.dMdecayingboxesnext_dM(Mnow,emsnow,alpha,2) Fun.dMdecayingboxesnext_dM(Mnow,emsnow,alpha,3) Fun.dMdecayingboxesnext_dM(Mnow,emsnow,alpha,4)];
        Fun.dMnext_dems = @(Mnow,emsnow,alpha) [ Fun.dMpermnext_dems(Mnow,emsnow,alpha) Fun.dMdecayingboxesnext_dems(Mnow,emsnow,alpha,2) Fun.dMdecayingboxesnext_dems(Mnow,emsnow,alpha,3) Fun.dMdecayingboxesnext_dems(Mnow,emsnow,alpha,4)];
        Fun.dMnext_dalpha = @(Mnow,emsnow,alpha) [ Fun.dMpermnext_dalpha(Mnow,emsnow,alpha) Fun.dMdecayingboxesnext_dalpha(Mnow,emsnow,alpha,2) Fun.dMdecayingboxesnext_dalpha(Mnow,emsnow,alpha,3) Fun.dMdecayingboxesnext_dalpha(Mnow,emsnow,alpha,4)];
        
    otherwise
        
        temp_matrix = Params.Mmatrix^Params.timestep;
        
        % index_nextsink identifies which sink are differentiating; output
        % is what are differentiating with repect to
        Fun.dMnext_dM = @(Mnow,emsnow,index_nextsink) bsxfun(@times,ones(size(Mnow,1),size(Mnow,2)),temp_matrix(index_nextsink,:)); % TxSinks; sinks index_nextsink is determined by multiplying row index_nextsink of temp_matrix by column vector of current reservoir values
        Fun.dMnext_dems = @(Mnow,emsnow,index_nextsink) ones(size(emsnow,1),1)*Params.emsinks(index_nextsink); % Tx1
end
switch Params.carbonmodel
    case 'dice' % atmospheric CO2 is the first reservoir
        Fun.dforcing_dM = @(M) [ Params.f2x*(1./M(:,1))/log(2)   zeros(size(M,1),size(M,2)-1) ]; % T x sinks
    otherwise % atmospheric CO2 is sum of all reservoirs
        Fun.dforcing_dM = @(M) (Params.f2x*(1./sum(M,2))/log(2))*ones(1,size(M,2)); % T x sinks
end

Fun.dTnext_dT = @(Tnow,Toceannow,Mnext) ones(size(Tnow,1),1) + ones(size(Tnow,1),1)*(Params.timestep/5)*Params.phi1*( -Params.lambda - Params.phi3 ); % Tx1
Fun.dTnext_dTocean = @(Tnow,Toceannow,Mnext) ones(size(Toceannow,1),1)*(Params.timestep/5)*Params.phi1*( Params.phi3 ); % Tx1
Fun.dTnext_dM = @(Tnow,Toceannow,Mnext) (Params.timestep/5)*Params.phi1*( Fun.dforcing_dM(Mnext) ); % T x sinks

Fun.dToceannext_dT = @(Tnow,Toceannow) ones(size(Tnow,1),1)*(Params.timestep/5)*Params.phi4; % Tx1
Fun.dToceannext_dTocean = @(Tnow,Toceannow) ones(size(Toceannow,1),1) - ones(size(Toceannow,1),1)*(Params.timestep/5)*Params.phi4; % Tx1

% Functions related to FAIR's carbon cycle feedbacks, from Dietz et al. (2020)
Fun.alpha_analytic = @(IRF) 0.0107*exp(0.0866*IRF); % analytic approximation to alpha; see Dietz et al. (2020) pg 37 equation (13)
Fun.IRF1 = @(T,M,cumulems) Params.R0 + Params.RT*T + Params.RC*(cumulems-(sum(M,2)-588)); % 100-year integrated impulse response target
Fun.IRF2 = @(alpha) bsxfun(@plus,100*Params.emsinks(1), sum( bsxfun(@rdivide,Params.emsinks(2:4),bsxfun(@rdivide,Params.decay(2:4),alpha)).*(1-exp(-100*bsxfun(@rdivide,Params.decay(2:4),alpha))) , 2) ); % Tx1, where alpha is Tx1 and Params.emsinks and Params.decay are 1x4; 100-year integrated impulse response implied by carbon model, for given alpha

Fun.dIRF1_dT = @(T,M,cumulems) Params.RT*ones(size(T,1),1); % Tx1
Fun.dIRF1_dM = @(T,M,cumulems) -Params.RC*ones(size(M,1),size(M,2)); % Tx4
Fun.dIRF1_dems = @(T,M,cumulems) Params.RC*ones(size(cumulems,1),1); % Tx1
Fun.dIRF2_dalpha = @(alpha) sum( bsxfun(@rdivide,Params.emsinks(2:4),bsxfun(@rdivide,Params.decay(2:4),ones(size(alpha,1),1))).*(1-exp(-100*bsxfun(@rdivide,Params.decay(2:4),alpha))) , 2)  ...
    + sum( bsxfun(@rdivide,Params.emsinks(2:4),bsxfun(@rdivide,Params.decay(2:4),alpha)).*(-100*bsxfun(@rdivide,Params.decay(2:4),alpha.^2)).*(exp(-100*bsxfun(@rdivide,Params.decay(2:4),alpha))) , 2) ; % Tx1


%% Unit conversions

Params.gtc_per_ppm = 2.16; % unit conversion
Params.co2_per_c = 44/12; % unit conversion


%% Abatement and marginal abatement cost, for calculating tax from abatement policy

Fun.abatement = @(sigma,abaterate,Ygross) sigma.*abaterate.*Ygross;
Fun.mac = @(sigma,abaterate,Ygross,t) Params.abate_exp*Fun.psi(t,sigma).*(abaterate.^Params.abate_exp).*Ygross./Fun.abatement(sigma,abaterate,Ygross) ...
    *(1e12/1e9)/Params.co2_per_c; % converted to 2010$/tCO2 from trillion 2010$ / Gt C


%% Exogenous trajectories
          
Params.discfactor = 1./((1+Params.discrate).^(Params.timestep*(0:Params.horizon-1)'));

Params.pop(1,1) = Params.L0;
Params.tfp(1,1) = Params.tfp0;
Params.sigma(1,1) = Params.sigma0;
Params.psi(1,1) = Fun.psi(1,Params.sigma(1,1));

for t=2:Params.horizon
    Params.pop(t,1) = Fun.popnext(Params.pop(t-1,1));
    Params.tfp(t,1) = Fun.tfpnext(t-1,Params.tfp(t-1,1));
    Params.sigma(t,1) = Fun.sigmanext(t-1,Params.sigma(t-1,1));
    Params.psi(t,1) = Fun.psi(t,Params.sigma(t,1));
end


