% This script runs the unit test for the computation of BN matrix

addpath('lib')
addpath('lib/util')
addpath('lib/propagator')
addpath('lib/measurements')
addpath('lib/inputs')

clear
clc
format long g

%% Set up parameters
loadDebugHW1;

% Time

dt = 60;
tf = 4*3600;
t  = t0+(0:dt:tf); % time vector to save s/c states in sec

numEpochs = length(t); 

gammas = zeros(numEpochs,1);
betas  = zeros(numEpochs,1);
alphas = zeros(numEpochs,1);

for ii = 1:numEpochs
        
    % Compute the rotation matrix between the inertial and body frame
    Wt  = wrapTo2Pi(Wt0 + Wdot*(t(ii)-t0)); % position of prime meridian in rad
    BN  = BodyFrame(Wt,DEC,RA);             % inertial to body frame rotation matrix
    NB  = BN';                              % body frame to inertial frame rotation matrix
 
    % Obtain 3-1-3 Euler angles
    [gamma,beta,alpha] = invBN(BN);
    
    gammas(ii) = gamma;
    betas(ii)  = beta;
    alphas(ii) = alpha;
    
end % For ii

figure(1);
plot(alphas*180/pi)
title('alpha angle')

figure(2);
plot(betas*180/pi)
title('beta angle')

figure(3);
plot(gammas*180/pi)
title('gamma angle')
