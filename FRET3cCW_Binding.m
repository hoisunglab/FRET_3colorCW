assert(false,'This script is not for the direct excution.');

% This script includes example of 3-color CW binding analysis.

%% Input data

% Example data.
p = pwd; %put your installation directory

input_files{1}=[p '\Binding\trj_DA12.mat'];
input_files{3}=[p '\Binding\trj_DA2.mat'];
input_files{4}=[p '\Binding\index_DA12.mat'];
input_files{6}=[p '\Binding\index_DA2.mat'];

input_data=input_bundle(input_files,'Binding');

burstbint3rall=input_data{1};
cumindex=input_data{2};
indexone=input_data{3};

%% pre-defined parameters 
chans=[1 3 4]; %[accpetor2 acceptor2 donor]
fparam={[0.05 0.0387]', [], [], 0.2}; %leakages and FRET efficiencies for acceptor dark states
                                      % fparam{1}=[donor->acceptor1
                                      % leakage, donor->acceptor2 leakage];
                                      % fparam{4}=acceptor1->acceptor2
                                      % leakage
                                      % fparam{2}=E1 FRET values for acceptor2 dark
                                      % states. When the model has n states,
                                      % you need n values.
                                      % fparam{2}=E2 FRET values for acceptor1
                                      % dark states. When the model has n
                                      % states, you need n values.
                                      % fparam{2}, fparam{3} will be
                                      % assigned later depending on models
                                      
isfastbinding=true;         % Analysis on binding data

%% Parameter initialization: 2-state analysis with acceptor blinkings
  
ubid=2;stateidML=[1 3]; % Kinetics model is State1-State2-State3-State1-... (all connected). State1: bound state (stateidML=1), State2: unbound state (ubid=2), State3: bound state (acceptor dark, stateidML=3)
isA2bleach=2;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model. Set this value to 2 when detailed balance is broken.
fparam{2}=[0 0]';
fparam{3}=[0 0]';
initparams=[0.06 0.35 0.7 0.1 0.65 0.55 0.02 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
                ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
                0.4 0.2 0.3 ...   % Labeling efficiency of A2, apparent labeling efficiency of A2, and rates for bleaching or blinking relaxation of A2 (kb + kd)
                0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
                0.35 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                20 20 0.9 0.9]'; % ks to bright state of A1 and A2, populations of the birght state of A1 and A2
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.4 0.9; 0.45 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.05 0.3; 0.01 1; ...
            0.1 5; ...
            0.1 0.6; ...
            5 200; 5 200; 0.8 0.999; 0.8 0.999];
        
% parameter optimization
resparams=mlhrateeiglingen3c1ExABIC(initparams,LUbounds,fparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
% error estimation. Covariance matrices are used to estimate errors.
drvintv=0.005; %relative perturbation for covariance matrices
drvpnt=5;  % number of perturbations on each side. Total 2*drvpnt perturbations will be evaluated around the determined parameters.

nstate=(length(initparams)-2-isA2bleach)/6; %total number of physical states
pfactor=cumprod([1; 1-resparams(end-nstate-2:end-5)]); %conversion factor for relative population

% Because there are two state for a bound state (stateML == [1 3]), convert
% k value for the bound state
% Error is calculated for k (sum of rates between adjacent states) and relative population p.
resparam=[resparams(1:5*nstate-1+isA2bleach); resparams(end-nstate-2:end-4).*pfactor; resparams(end-3:end)]; % Convert f to population p
[errorparams, logmlh, bic]=mlhrateeiglingen3c1ExABIerrorC(resparam,fparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans,drvintv,drvpnt);
%% Example1, without acceptor blinkings. DA1A2 segments only
 
% n-state linear model
% 2 state
ubid=2;stateidML=1; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
initparams=[0.06 0.35 0.7 0.05 0.65 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.35]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
               % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.5 0.8; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; ...
            0.1 0.5];
        
%     isA2bleach=0;    % False: A2 bright and dark state in the bound state is not connected. 
%                          % So, the transition occur through dissociation and rebinding. This will increase associateion/dissociation rates.
%     initparams=[0.06 0.35 0.7 0.05 0.65 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
%                 ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
%                 0.4 ...   % Labeling efficiency of A2
%                 0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
%                 0.35]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
%                    % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
%     LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.5 0.8; ...
%                 0.2 0.8; ...
%                 0.1 5; ...
%                 0.1 0.5];

% 3 off state (B - U - B')
ubid=2;stateidML=1; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
initparams=[0.06 0.35 0.2 ...   % (1:n): E1app(1) ~ E1app(n)
            0.7 0.05 0.45 ...   % (n+1:2n): E2app(1) ~ E2app(n)
            0.7 0.5 ...   % (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound)
            ...             % If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...             % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.33 0.8]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
               % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.01 0.2; 0.2 0.5; 0.1 0.4; ...
            0.4 0.9; 0.01 0.2; 0.2 0.6; ...
            0.5 0.8; 0.4 0.7; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; 0.2 5;
            0.1 0.5; 0.6 0.95];

% 4 state (B - I - U - B')
ubid=3;stateidML=1; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
initparams=[0.06 0.2 0.35 0.2 ...   % (1:n): E1app(1) ~ E1app(n)
            0.7 0.45 0.05 0.45 ...   % (n+1:2n): E2app(1) ~ E2app(n)
            0.7 0.5 0.5 ...   % (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound)
            ...             % If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...             % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            10 0.2 2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.33 0.1 0.8]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
               % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.01 0.2; 0.1 0.4; 0.2 0.5; 0.1 0.4; ...
            0.4 0.9; 0.2 0.6; 0.01 0.2; 0.2 0.6; ...
            0.5 0.8; 0.4 0.7; 0.4 0.7; ...
            0.2 0.8; 0.001 0.2; ...
            1 50; 0.1 5; 0.2 5;
            0.1 0.5; 0.01 0.3; 0.6 0.95];

%for the above models only DA1A2 segments (burstbint3rall(1)) are used.
resparams=mlhrateeiglingen3c1ExC(initparams,LUbounds,burstbint3rall(1),cumindex(1),indexone(1),isfastbinding,isA2bleach,ubid,stateidML,chans);

%% Example2. Without acceptor blinkings. All segments including acceptor bleached segments
ubid=2;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
initparams=[0.06 0.35 0.7 0.1 0.65 0.55 0.02 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.35]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
               % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.4 0.9; 0.45 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; ...
            0.1 0.6];

ubid=2;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=2;    % 1: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
                % 2: no detailed balance.
initparams=[0.075 0.34 0.72 0.1 0.7 0.57 0.02 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.2 0.3 ...   % Labeling efficiency of A2, apparent labeling efficiency of A2, and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.35]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
               % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.4 0.9; 0.45 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.05 0.3; 0.01 1; ...
            0.1 5; ...
            0.1 0.6];

% 3 linear state (B - I - U)
ubid=3;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=2;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
initparams=[0.07 0.07 0.35 ...   % (1:n): E1app(1) ~ E1app(n)
            0.8 0.6 0.1 ...   % (n+1:2n): E2app(1) ~ E2app(n)
            0.75 0.5 ...   % (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound)
            0.6 0.4 0.03 ...             % If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...             % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.3 0.2 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            1 1 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.33 0.3 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
            ]';
LUbounds=[0.01 0.2; 0.03 0.4; 0.2 0.5; ...
            0.4 0.9; 0.2 0.75; 0.01 0.2; ...
            0.5 0.8; 0.4 0.7; ...
            0.4 0.8; 0.15 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.05 0.5; 0.01 1; ...
            0.1 10; 0.1 10;
            0.1 0.7; 0.1 0.7];

% 3 off state (B - U - B')
ubid=2;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
initparams=[0.06 0.35 0.2 ...   % (1:n): E1app(1) ~ E1app(n)
            0.7 0.05 0.45 ...   % (n+1:2n): E2app(1) ~ E2app(n)
            0.7 0.5 ...   % (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound)
            0.55 0.02 0.3 ...             % If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...             % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.33 0.8]';    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
               % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.01 0.2; 0.2 0.5; 0.1 0.4; ...
            0.4 0.9; 0.01 0.2; 0.2 0.6; ...
            0.5 0.8; 0.4 0.7; ...
            0.4 0.8; 1e-6 0.1; 0.15 0.5; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; 0.2 5;
            0.1 0.5; 0.6 0.95];
            
resparams=mlhrateeiglingen3c1ExC(initparams,LUbounds,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
    
%% parameter conversion and error estimation for example 1 and 2
if isfastbinding,
    if all (stateidML == [1 3]),
        nstate=(length(initparams)+2-isA2bleach)/6;
    else
        nstate=(length(initparams)+2-isA2bleach)/5;
    end
else
    nstate=(length(initparams)+2)/(length(stateidML)+3);
end

pfactor=cumprod([1; 1-resparams(end-nstate+2:end-1)]);

drvintv=0.005; %relative perturbation for covariance matrices
drvpnt=5;  % number of perturbations on each side. Total 2*drvpnt perturbations will be evaluated around the determined parameters.

% resparam conversion for stateML == 1 (example 1)
% Error is calculated for k (sum of rates between adjacent states) and relative population p.
resparam=[resparams(1:4*nstate-1+isA2bleach); resparams(end-nstate+2:end).*pfactor]; % Convert f to population p
[errorparams, logmlh, bic]=mlhrateeiglingen3c1ExerrorC(resparam, burstbint3rall(1),cumindex(1),indexone(1),isfastbinding,isA2bleach,ubid,stateidML,chans,drvintv,drvpnt);

% resparam conversion for stateML == [1 3] (example 2)
% Error is calculated for k (sum of rates between adjacent states) and relative population p.
resparam=[resparams(1:5*nstate-1+isA2bleach); resparams(end-nstate+2:end).*pfactor]; % Convert f to population p, include DA2
[errorparams, logmlh, bic]=mlhrateeiglingen3c1ExerrorC(resparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans,drvintv,drvpnt);

%% Example 3, with acceptor blinkings. DA1A2 only.
    
% n-state linear model with acceptor blinking
% 2 state
ubid=2;stateidML=1; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
fparam={[0.05 0.0387]', [0 0]', [0.6 0.0387]', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
                                % [effb1 effb2 effb12] are donor and accetor leak into other channels (effb1 = A1/(A1+D), effb2 = A2/(A2+D), effb12 = A2/(A1+A2)
                                % effs2 and effs3 are FRET efficiencies in state 2 (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A2)). These values are fitting parameters in some conditions
initparams=[0.06 0.35 0.7 0.05 0.65 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.35 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
            20 20 0.9 0.9]'; % ks to bright state of A1 and A2, populations of the birght state of A1 and A2
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.5 0.8; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; ...
            0.1 0.5; ...
            5 100; 5 100; 0.8 0.999; 0.8 0.999];

resparams=mlhrateeiglingen3c1ExABIC(initparams,LUbounds,fparam,burstbint3rall(1),cumindex(1),indexone(1),isfastbinding,isA2bleach,ubid,stateidML,chans);

%% Example 4, with acceptor blinkings. All segments including acceptor bleached segments
% 2 state
ubid=2;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
fparam={[0.05 0.0387]', [0 0]', [0 0]', 0.2};
initparams=[0.06 0.35 0.7 0.1 0.65 0.55 0.02 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.35 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
            20 20 0.9 0.9]'; % ks to bright state of A1 and A2, populations of the birght state of A1 and A2
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.4 0.9; 0.45 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; ...
            0.1 0.6; ...
            5 100; 5 100; 0.8 0.999; 0.8 0.999];

% 2 state
ubid=2;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=2;    % 1: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
                % 2: no detailed balance.
fparam={[0.05 0.0387]', [0 0]', [0 0]', 0.2};
initparams=[0.06 0.35 0.7 0.1 0.65 0.55 0.02 ... % (1:n): E1app(1) ~ E1app(n); (n+1:2n): E2app(1) ~ E2app(n); (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound). If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...                        % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.2 0.3 ...   % Labeling efficiency of A2, apparent labeling efficiency of A2, and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.35 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
            20 20 0.9 0.9]'; % ks to bright state of A1 and A2, populations of the birght state of A1 and A2
LUbounds=[0.01 0.2; 0.2 0.5; 0.4 0.9; 0.01 0.2; 0.4 0.9; 0.45 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.05 0.3; 0.01 1; ...
            0.1 5; ...
            0.1 0.6; ...
            5 200; 5 200; 0.8 0.999; 0.8 0.999];

% 3 off state (B - U - B')
ubid=2;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
fparam={[0.05 0.0387]', [0 0 0]', [0 0 0]', 0.2};
initparams=[0.06 0.35 0.2 ...   % (1:n): E1app(1) ~ E1app(n)
            0.7 0.05 0.45 ...   % (n+1:2n): E2app(1) ~ E2app(n)
            0.7 0.5 ...   % (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound)
            0.55 0.02 0.3 ...             % If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...             % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.05 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            0.2 2 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.33 0.8 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
            20 20 0.9 0.9]'; % ks to bright state of A1 and A2, populations of the birght state of A1 and A2
LUbounds=[0.01 0.2; 0.2 0.5; 0.1 0.4; ...
            0.4 0.9; 0.01 0.2; 0.2 0.6; ...
            0.5 0.8; 0.4 0.7; ...
            0.4 0.8; 1e-6 0.1; 0.15 0.5; ...
            0.2 0.8; 0.001 0.2; ...
            0.1 5; 0.2 5;
            0.1 0.5; 0.6 0.95; ...
            5 100; 5 100; 0.8 0.999; 0.8 0.999];

% 3 linear state (B - I - U)
ubid=3;stateidML=[1 3]; %  For isfastbinding == true
isA2bleach=2;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
fparam={[0.05 0.0301]', [0 0 0]', [0 0 0]', 0.2};
initparams=[0.07 0.07 0.35 ...   % (1:n): E1app(1) ~ E1app(n)
            0.8 0.6 0.1 ...   % (n+1:2n): E2app(1) ~ E2app(n)
            0.7 0.5 ...   % (2n+1:3n-1): E1(1) (= (A1+A2)/n) ~ E1(n) except E1(unbound), which is the same as E1app(unbound)+E2app(unbound)
            0.7 0.5 0.03 ...             % If stateidML == [1 3], (3n:4n-1): E2(1) - E2(n)
            ...             % For isfastbinding == false, (2n+1:4n) can be E1 or E2 or both depending on the stateidML 1 or [1 2] or [1 3] or [1 2 3].
            0.4 0.3 0.2 ...   % Labeling efficiency of A2 and rates for bleaching or blinking relaxation of A2 (kb + kd)
            1 1 ...   % k(1) ~ k(n-1), k(1) is the sum of foward and backward rates between state 1 and 2;
            0.33 0.3 ...    % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
            20 20 0.95 0.9]'; % ks to bright state of A1 and A2, populations of the birght state of A1 and A2
LUbounds=[0.01 0.2; 0.03 0.4; 0.2 0.5; ...
            0.4 0.9; 0.2 0.8; 0.01 0.2; ...
            0.5 0.85; 0.4 0.7; ...
            0.4 0.8; 0.15 0.7; 1e-6 0.1; ...
            0.2 0.8; 0.05 0.5; 0.01 1; ...
            0.1 10; 0.1 10;
            0.1 0.7; 0.1 0.7; ...
            5 200; 5 200; 0.8 0.999; 0.8 0.999];
            
resparams=mlhrateeiglingen3c1ExABIC(initparams,LUbounds,fparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
            
%% parameter conversion and error estimation for example 3 and 4
if isfastbinding,
    if all (stateidML == [1 3]),
        nstate=(length(initparams)-2-isA2bleach)/6;
    else
        nstate=(length(initparams)-2-isA2bleach)/5;
    end
else
    nstate=(length(initparams)-2)/(length(stateidML)+3);
end

pfactor=cumprod([1; 1-resparams(end-nstate-2:end-5)]);

drvintv=0.005; %relative perturbation for covariance matrices
drvpnt=5;  % number of perturbations on each side. Total 2*drvpnt perturbations will be evaluated around the determined parameters.

% resparam conversion for stateML == 1 (Example 3)
% Error is calculated for k (sum of rates between adjacent states) and relative population p.
resparam=[resparams(1:4*nstate-1+isA2bleach); resparams(end-nstate-2:end-4).*pfactor; resparams(end-3:end)]; % Convert f to population p
[errorparams, logmlh, bic]=mlhrateeiglingen3c1ExABIerrorC(resparam,fparam,burstbint3rall(1),cumindex(1),indexone(1),isfastbinding,isA2bleach,ubid,stateidML,chans,drvintv,drvpnt);

% resparam conversion for stateML == [1 3] (Example 4)
% Error is calculated for k (sum of rates between adjacent states) and relative population p.
resparam=[resparams(1:5*nstate-1+isA2bleach); resparams(end-nstate-2:end-4).*pfactor; resparams(end-3:end)]; % Convert f to population p
[errorparams, logmlh, bic]=mlhrateeiglingen3c1ExABIerrorC(resparam,fparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans,drvintv,drvpnt);