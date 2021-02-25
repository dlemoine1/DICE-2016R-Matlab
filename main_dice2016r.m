% Replicates DICE-2016R

% Last edited: February 24, 2021 by Derek Lemoine

clear all;
close all;

global iteration_number


%% User Options %%%%%%%%%%%%%%%

customdir = '2020-09'; % prefix for save directory; not used if running on HPC

Params.timestep = 5; % timestep, in years; default: 5
Params.numyears = 500; % number of years in horizon; default: 500

Params.climatemodel = 'dice'; % 'dice' or 'update'; determines whether use DICE-2016R climate dynamics or Geoffroy et al dynamics (see Dietz et al. 2020)
Params.carbonmodel = 'dice'; % 'dice','gol','joos','fair'; determines whether use DICE-2016R carbon dynamics or one of the other models (see Dietz et al. 2020)
Params.damagemodel = 'dice'; % which damage distribution: 'dice' is DICE-2016R, and 'expert' is Pindyck (2019) damages as implemented by Lemoine (2021)

Params.fixsavings = -99; % <0 (default): endogenize savings rate; >=0: savings rate fixed at this value
Params.dicenegems = 1; % =1 (default): use DICE's constraint allowing negative emissions; =0: do not allow negative emissions
Params.dicecumulems = 1; % =1 (default): use DICE's constraint on cumulative emissions; =0: no constraint
Params.dicelrsavings = 1; % =1 (default): fix long-run savings rate as in DICE; =0: don't
Params.dicefirstperabate = 0; % =1: fix first-period abatement as in DICE, so optimization begins only from second-period; =0 (default): optimize all periods
Params.optimizeonlysavings = 0; % =1: optimize only savings and then calculate the social cost of carbon; =0 (default): optimize abatement as well

% Computational options
Params.transitionsasconstraints = 1; % =0: solve by guessing policy and simulating trajectories; =1 (default): solve by also guessing states and treating transition equations as constraints
Params.consumptionascontrol = 0; % =1: consumption is a control and output is a constraint; =0 (default): consumption is deduced from output; only relevant if Params.transitionsasconstraints=1
Params.scaleconstraints = 1; % =1 (default): scale the constraints to be fractional deviations; =0: don't
Params.scalevars = 1; % =1 (default): scale the variables by the initial guess; =0: don't
Params.dohpc = 0; % =0: Are running on personal computer; =1: Are running remotely on high performance computing (HPC)
Params.useknitro = 0; % =1: use Knitro; =0: use fmincon
Params.knitroalgorithm = 3; % algorithm to use with Knitro, can use 1-4; default: 3
Params.fminconalgorithm = 'interior-point'; % only relevant if Params.useknitro=0 or Params.dohpc=1; default: 'interior-point'
Params.parallelize = 0; % Only relevant when Params.dohpc==1 and Params.transitionsasconstraints~=1.  =1: Parallelize gradients in fmincon.  =0: Don't
Params.screenreport = 0; % =1: report output summary to screen as optimize; =0: don't

% Obtain directory
[maindir,~,~] = fileparts(mfilename('fullpath'));
addpath(maindir);

%Exporting Command Window Outputs to .txt file
cd(maindir);
diary 'FileRunOutput.txt'



%% Parameterization %%%%%%%%%%%%%%%

% Run parameterization script
run sub_parameters

% More computational parameters
if Params.dohpc==1 % HPC doesn't have knitro license so set up to do fmincon
    Params.useknitro = 0;
end
if Params.transitionsasconstraints==1
    if Params.useknitro
        Params.knitro_options_file = ['dice2016r_knitro_alg' num2str(Params.knitroalgorithm) '.opt'];
        Params.optimize_options_list = [];
    else
        Params.optimize_options_list = optimoptions(@fmincon,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',1e4,'MaxIterations',5e3,'Display','final-detailed','Algorithm',Params.fminconalgorithm,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
        %Params.optimize_options_list = optimset('GradObj','on','GradConstr','on','DerivativeCheck','on');    
    end    
else
    Params.knitro_options_file = [];
    if Params.dohpc==1 && Params.parallelize==1 % Parallelize gradient calculations
        parpool(5);
        Params.optimize_options_list = optimoptions(@fmincon,'Display','off','Algorithm',Params.fminconalgorithm,'UseParallel',true);
    else
        Params.optimize_options_list = optimoptions(@fmincon,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',3e4,'MaxIterations',5e3,'Display','final-detailed','Algorithm',Params.fminconalgorithm);
    end
end




%% Error check %%%%

if Params.fixsavings>=1
    error('Savings rate too large');
elseif Params.fixsavings>=0
    disp(['Fixing savings rate at ' num2str(Params.fixsavings)]);
end
if Params.dicelrsavings~=1 && Params.dicelrsavings~=0
    error('Params.dicelrsavings not binary.');
end

switch Params.climatemodel
    case 'dice'
        disp('DICE-2016R climate dynamics');
    case 'update'
        disp('Updated climate dynamics');
    otherwise
        error('Unrecognized climate model')
end

switch Params.carbonmodel
    case 'dice'
        disp('DICE-2016R carbon dynamics');
    case 'gol'
        disp('Golosov et al carbon dynamics');
    case 'joos'
        disp('Joos et al carbon dynamics');
    case 'fair'
        disp('FAIR carbon dynamics, including endogenous feedback parameter');
    otherwise
        error('Unrecognized carbon model')
end

switch Params.damagemodel
    case 'dice'
        disp('Damages from DICE-2016R');
    case 'expert'
        disp('Damages based on Pindyck (2019) expert survey, via Lemoine (2021)');
    otherwise
        error('Unrecognized damage type');
end


%% Directory Structure %%%%%%%%%%%%%%%

if Params.dohpc ~= 1 % set up directory for a run on personal computer
    
    Params.savedir = [maindir filesep 'output' filesep customdir ];    
    
    Params.savedir = [Params.savedir '_climate-' Params.climatemodel '_carbon-' Params.carbonmodel '_dam-' Params.damagemodel];
    Params.filenaming = ['cli-' Params.climatemodel '_car-' Params.carbonmodel '_dam-' Params.damagemodel];
    
    if Params.fixsavings>=0
        Params.savedir = [Params.savedir '_fixedsavings-' num2str(Params.fixsavings)];
        Params.filenaming = [Params.filenaming '_fixedsavings-' num2str(Params.fixsavings)];
    end        
    
    if Params.dicenegems~=1
        Params.savedir = [Params.savedir '_nonegems'];
        Params.filenaming = [Params.filenaming '_nonegems'];
    end
    
    if Params.dicecumulems~=1
        Params.savedir = [Params.savedir '_nocumulems'];
        Params.filenaming = [Params.filenaming '_nocumulems'];
    end
    
    if Params.dicelrsavings~=1
        Params.savedir = [Params.savedir '_freelrsavings'];
        Params.filenaming = [Params.filenaming '_freelrsavings'];
    end
       
    if Params.optimizeonlysavings==1
        Params.savedir = [Params.savedir '_noabate'];
        Params.filenaming = [Params.filenaming '_noabate'];
    end
    if Params.dicefirstperabate==1
        Params.savedir = [Params.savedir '_exogabate1'];
        Params.filenaming = [Params.filenaming '_exogabate1'];
    end    
    
    if Params.transitionsasconstraints==1
        Params.savedir = [Params.savedir '_constraint-sol'];
    else
        Params.savedir = [Params.savedir '_trajectory-sol'];
    end
    
    if Params.useknitro==1
        Params.savedir = [Params.savedir '_knitro'];
    else
        Params.savedir = [Params.savedir '_fmincon'];
    end
    
    Params.savedir = [Params.savedir filesep];
    
    if ~exist(Params.savedir, 'dir')
        mkdir(Params.savedir);
    end
    cd(Params.savedir);
    
    if Params.useknitro==1 && ~isempty(Params.knitro_options_file)
        copyfile([maindir filesep Params.knitro_options_file],Params.savedir);
    end
    
else % don't have permission to mkdir from Matlab on HPC, so use same directory for saving output
    
    Params.savedir = maindir;
    
end

disp(['Saving to ' Params.savedir]);

if Params.dohpc==1
    delete([Params.savedir 'iteration_report.txt']);
end



%% Optimization %%%%%%%%%%%%%%%

% load guess and variable bounds (ub and lb)
run sub_loadguesses;

% initialize iteration tracking
iteration_number = 0;

% objective and nonlinear constraints
objective = @(x) utilityobjective(x, 1, Params.horizon, Fun, Params);
nonlcon = @(x) nonlcon_utilmax(x, 1, Params.horizon, Fun, Params);
if Params.fixsavings<0
    disp('Beginning to optimize abatement and savings rates.');
else
    disp('Beginning to optimize abatement rate.');
end


% optimize
dosolve=1; % when = x, will try at most x times
elapsedtime = 0;
tic;
while dosolve>0
    if Params.useknitro == 1
        [out_controls,fval,exitflag,outputstruct,lambdastruct,gradient,hessian] = knitromatlab(objective,guess(:),[],[],[],[],lb(:),ub(:),nonlcon,[],Params.optimize_options_list,Params.knitro_options_file);
        if exitflag<0
            disp(['Warning: Optimization may have failed, with exitflag ' num2str(exitflag)]);
            dosolve = dosolve - 1;
            if dosolve > 0
                guess = out_controls;
            end
        else
            disp('Solved successfully')
            dosolve = 0;
        end
    else
        [out_controls,fval,exitflag,outputstruct,lambdastruct,gradient,hessian] = fmincon(objective,guess(:),[],[],[],[],lb(:),ub(:),nonlcon,Params.optimize_options_list);
        if exitflag<=0
            disp(['Warning: Optimization may have failed, with exitflag ' num2str(exitflag)]);
            dosolve = dosolve - 1;
            if dosolve > 0
                guess = out_controls;
                % try sqp algorithm
                Params.optimize_options_list = optimoptions(@fmincon,'Algorithm','sqp','StepTolerance',1e-10,'ConstraintTolerance',1e-12,'MaxFunctionEvaluations',1e4,'MaxIterations',5e3,'Display','final-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
            end
        else
            disp('Solved successfully')
            dosolve = 0;
        end
    end
    elapsedtime = elapsedtime + toc;
    disp(['Time spent optimizing: ' num2str(elapsedtime/60) ' minutes.']);
end

% store constraints
[out_cneq,out_ceq] = nonlcon(out_controls);

if Params.transitionsasconstraints
   
    disp(['Max abs equality constraint is ' num2str(max(abs(out_ceq)))]);
    
    % turn back into matrix
    out_controls = reshape(out_controls,Params.horizon,[]);
    
    % make sure state variable bounds didn't bind
    test_ub = ub;
    test_lb = lb;
    test_ub(:,[Params.col_abaterate Params.col_savingsrate]) = Inf; % don't care about controls binding, so get rid of that
    test_lb(:,[Params.col_abaterate Params.col_savingsrate]) = -Inf; % don't care about controls binding, so get rid of that
    ub_binding = sum( out_controls>=test_ub | abs(out_controls-test_ub)<=1e-4 , 1);
    lb_binding = sum( out_controls<=test_lb | abs(out_controls-test_lb)<=1e-4 , 1);
    if sum(ub_binding) > 0
        disp(['Upper bounds on states bind ' mat2str(ub_binding) ' times']);
    end
    if sum(lb_binding) > 0
        disp(['Lower bounds on states bind ' mat2str(lb_binding) ' times']);
    end
    clear test_ub test_lb;
            
end


%% Take Results %%%%%%%%%%%%%%%

% store the policy variables
if Params.transitionsasconstraints~=1
    
    % undo normalization
    out_controls = out_controls(:).*Params.normalization(:);
    
    out_policy = out_controls(:);
    
    out_policy2 = out_policy;
    if strcmp(Params.carbonmodel,'fair')
        alpha = out_policy2(1:Params.horizon,1);
        out_policy2(1:Params.horizon,:) = [];
    end
    abaterate = out_policy2(1:Params.horizon,1); % fraction
    if Params.fixsavings<0
        savingsrate = out_policy2(Params.horizon+1:end,1);
        if length(savingsrate)<Params.horizon
            savingsrate(end+1:Params.horizon,1) = Params.lrsavingsrate;
        end
    else
        savingsrate = Params.fixsavings*ones(Params.horizon,1);
    end
    clear out_policy2;
    
else % turn into appropriate vector
    
    % undo normalization
    out_controls = out_controls.*Params.normalization;
    
    if Params.fixsavings>=0
        if strcmp(Params.carbonmodel,'fair')
            out_policy = [ out_controls(:,Params.col_alpha); out_controls(:,Params.col_abaterate) ];
        else
            out_policy = out_controls(:,Params.col_abaterate);
        end
    else
        if strcmp(Params.carbonmodel,'fair')
            out_policy = [ out_controls(:,Params.col_alpha); out_controls(:,Params.col_abaterate); out_controls(:,Params.col_savingsrate) ];
        else
            out_policy = [ out_controls(:,Params.col_abaterate); out_controls(:,Params.col_savingsrate) ];
        end
        if Params.dicelrsavings==1
            out_policy(end-round(5*10/Params.timestep)+1:end,:) = [];
        end                    
    end
    
    abaterate = out_controls(:,Params.col_abaterate);
    if Params.fixsavings<0
        savingsrate = out_controls(:,Params.col_savingsrate);
        if length(savingsrate)<Params.horizon
            savingsrate(end+1:Params.horizon,1) = Params.lrsavingsrate;
        end
    else
        savingsrate = Params.fixsavings*ones(Params.horizon,1);
    end
    if strcmp(Params.carbonmodel,'fair')
        alpha = out_controls(:,Params.col_alpha);
    end    
end

[ C, Ynet, K, T, M, emsind, Tocean, Ygross, abatecost] = trajectory( out_policy, 1, Params.horizon, Fun, Params );

year = 2015 + [0:Params.timestep:Params.numyears-Params.timestep]';

Welfare = sum(Params.discfactor.*Params.pop.*Fun.utility(C,Params.pop));

abatement = Fun.abatement(Params.sigma,abaterate,Ygross); % Gt C

% cumulative industrial emissions
cumulemsind = cumsum(emsind(1:Params.horizon));
cumulemsind = Params.cumulems0 + [0; cumulemsind(1:end-1,:)];

% cumulative fossil use, in Gt C
cumulfossil = cumsum(max(0,emsind));
cumulfossil = Params.cumulems0 + [0; cumulfossil(1:end-1,:)];

% implied tax on emissions, from marginal abatement cost
emtax_pertCO2 = Fun.mac(Params.sigma,abaterate,Ygross,[1:Params.horizon]'); % 2010$/tCO2
% note that have to adjust last argument if not wanting to run from first period

% atmospheric CO2, in ppm
switch Params.carbonmodel
    case 'dice' % atmospheric CO2 is the first reservoir
        Carbon_ppm = M(:,1)/Params.gtc_per_ppm;
    otherwise % atmospheric CO2 is sum of all reservoirs
        Carbon_ppm = sum(M,2)/Params.gtc_per_ppm;
end

% approximate implied alpha
if strcmp(Params.carbonmodel,'fair')
    otherems = [100;Fun.otherems([1:Params.horizon-1]')]; % Params.horizon x 1 vector of cumulative emissions from deforestation (Gt C); is the leading 100 Gt C meant to adjust for all emissions prior to 2015?
    alphaimplied = Fun.alpha_analytic( Fun.IRF1(T,M,cumulemsind + cumsum(otherems)) );
end

% format for future guess
output_for_guess(:,Params.col_abaterate) = abaterate;
output_for_guess(:,Params.col_savingsrate) = savingsrate;
if Params.consumptionascontrol
    output_for_guess(:,Params.col_C) = C;
end
output_for_guess(:,Params.col_K) = K;
output_for_guess(:,Params.col_T) = T;
output_for_guess(:,Params.col_Tocean) = Tocean;
output_for_guess(:,Params.col_M) = M;
if strcmp(Params.carbonmodel,'fair')
   output_for_guess(:,Params.col_alpha) = alpha;
end

save(['workspace.mat']);


%% Calculate Social Cost of Carbon, based on code by both Lemoine and Healy %%%%%%%%%%%%%%%

%Predefine SCC array
SCC_pertCO2 = zeros(Params.horizon,1) ;

%Save initial welfare result
W0 = Welfare;

for tstep = 1:Params.horizon

    %Increment emissions in period tstep
    Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + (t==tstep);

    %Call trajectory function to calculate new consumption vector and welfare
    [ C_margems, ~, ~, ~, ~, ~, ~, ~, ~] = trajectory(out_policy, 1, Params.horizon, Fun, Params );
    
    %Calculate marginal welfare
    W1e = sum(Params.discfactor.*Params.pop.*Fun.utility(C_margems,Params.pop)); 
    
    %Calculate partial derivative of welfare with respect to consumption
    dWc = Params.discfactor(tstep)*Params.pop(tstep)*Fun.dutility_dC(C(tstep),Params.pop(tstep));
    
    %Calculate SCC
    SCC_pertCO2(tstep,1) = -((W1e - W0)/dWc)*1000/Params.co2_per_c ;  

end

%Restore otherems function
Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5);

clear W1e W0 dWc C_margems;


%% Output Results to Excel File, from Ben Healy %%%%%%%%%%%%%%%

%Run Output Script, from Ben Healy
run OutputResults


%% Final save

save(['workspace.mat']);
