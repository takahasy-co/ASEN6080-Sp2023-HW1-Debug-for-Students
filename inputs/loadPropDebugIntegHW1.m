% Integration tolerances
prop.opts     = odeset('abstol',1e-16,'reltol',3e-14);
prop.partials = 1;

% Time

dt = 60;
tf = 4*3600;
T  = t0+(0:dt:tf); % time vector to save s/c states in sec

% Bias params

prop.numBias   = 7;
prop.biasNames = {'GM','C20','C21','C22','C30','S21','S22'};
prop.biasIdx   = [7,8,9,10,11,12,13];
bias0          = [pars.Body.GM;
                   pars.Body.SH.C(3,1); % C20
                   pars.Body.SH.C(3,2); % C21
                   pars.Body.SH.C(3,3); % C22
                   pars.Body.SH.C(4,1); % C30
                   pars.Body.SH.S(3,2); % S21
                   pars.Body.SH.S(3,3); % S22
                 ];

% Dynamic params

prop.numDynamic   = 0;
prop.dynamicNames = {};
prop.dynamicIdx   = [];
dynamic0          = [
                    ];

% Stochastic params

prop.numStoch   = 0;
prop.stochNames = {};
prop.stochIdx   = [];
stoch0          = [
                  ];
                  
prop.numParams  = 6 + prop.numBias + prop.numDynamic + prop.numStoch;

% Initial state

x0 = [xsc0; dynamic0; bias0; stoch0; reshape(eye(prop.numParams),prop.numParams^2,1)];
