% This script verifies the STM computation by finite differencing.

addpath('lib')
addpath('lib/util')
addpath('lib/propagator')
addpath('lib/measurements')
addpath('inputs')

clear
clc
format long g

%% Set up parameters
loadDebugHW1;

% Propagator setup
loadPropDebugSTMFiniteDifferencingHW1;

% create trajectory
[t,y] = ode45(@(t,x) sc_EOM_Phi(t,x,pars,prop),T,x0,prop.opts);

plot3(y(:,1),y(:,2),y(:,3)); axis equal

y_nom_final = y(end,1:prop.numParams);
Phi_nominal = reshape(y(end,prop.numParams+1:end),prop.numParams,prop.numParams);

% Finite differencing.
% No partials necessary.

prop.partials = 0;

Phi_fd = zeros(prop.numParams,prop.numParams);

for ii = 1:prop.numParams
    
    fprintf('Param %d\n',ii)
    
    if ii <= 3
        
        delta = 1e-5;
        
    elseif ii <= 6
        
        delta = 1e-8;

    elseif ii == 7
        
        delta = 1e-3;

    else
        
        delta = 1e-8;
        
    end % For if
    
    x0_p     = x0;
    x0_p(ii) = x0_p(ii) + delta;
    
    [t_p,y_p] = ode45(@(t,x) sc_EOM_Phi(t,x,pars,prop),T,x0_p,prop.opts);
    
    y_nom_final  = y(end,1:prop.numParams);
    y_p_final    = y_p(end,1:prop.numParams);
    
    Phi_fd(:,ii) = (y_p_final - y_nom_final)'/delta;
       
end % For ii

Phi_nominal
Phi_fd

[deltaM, deltaMPct] = diffMatrix(Phi_nominal,Phi_fd);

deltaMPct

%% EOM
function dx = sc_EOM_Phi(t,x,pars,prop)

    % GM + SH
    dx = accel_hw1(t,x,pars.Body,prop);
    
end % For function
