% Body parameters
pars.Body.GM        = 3.98600435e5; % Gravitational parameter km^3/s^2
pars.Body.Rref      = 6378.136; % reference radius for spherical harmonics in km
pars.Body.t0        = 0; % epoch for prime meridian
pars.Body.Wt0       = 32.9*pi/180; % location of prime meridian at epoch in rad
pars.Body.Wdot      = 2*pi/(86164); % rotation rate of prime meridian in rad/s
pars.Body.h_hat     = pars.Body.Wdot*[0;0;1]; % angular velocity vector of body fixed frame relative to the inertial frame in rad/s 
pars.Body.RA        = 83.4*pi/180; % right ascension of Body pole in rad
pars.Body.DEC       = 88.5*pi/180; % declination of Body pole in rad
                    
pars.Body.SH.degree = 3; % The degree and order of the harmonics.

pars.Body.SH.C      = zeros(pars.Body.SH.degree+1,pars.Body.SH.degree+1);
pars.Body.SH.C(1,1) = 1.0;
pars.Body.SH.C(3,1) = -1.082635e-3;
pars.Body.SH.C(3,2) = -2.58198e-10;
pars.Body.SH.C(3,3) = 1.57457e-6;
pars.Body.SH.C(4,1) = 2.532473e-6;

pars.Body.SH.S      = zeros(pars.Body.SH.degree+1,pars.Body.SH.degree+1);
pars.Body.SH.S(3,2) = 1.80326e-9;
pars.Body.SH.S(3,3) = -9.03867e-7;

pars.Body.rho0      = 3.614e-4; % atmospheric density reference in kg/km^3
pars.Body.r0        = 700+pars.Body.Rref; % atmospheric reference radisu in km
pars.Body.H         = 88.667; % km
pars.Body.Stations.State = [
                              [-5127.51;
                              -3794.16;
                              0;
                              zeros(3,1)],...
                              [3860.91;
                              3238.49;
                              3898.09;
                              zeros(3,1)],...
                              [549.505;
                              -1380.872;
                              6182.197;
                              zeros(3,1)]
                              ]; % station state vectors on Body (body frame) in km and km/s
pars.Body.Stations.IDs = [101;337;394]; % station IDs

% s/c parameters
pars.SC.mass = 970;  % s/c mass in kg
pars.SC.Area = 3e-6; % s/c cross-sectional area in km^2
pars.SC.Cd   = 2;    % s/c coefficient of drag
pars.SC.CR   = 1;    % s/c coefficient of reflectivity

% Astro constants
pars.Astro.speedOfLight = 2.99792458e5;    % [km/s]        Speed of light
pars.Astro.AU           = 1.495978707e+08; % [km]          Mean distance between the Sun and the Earth
pars.Astro.boltzman     = 5.67e-8;         % [w/(m^2*k^4) = kg/(s^3*k^4)] Stefan-Boltzman's constant
pars.Astro.R_sun        = 6.96000e5;       % [km]          Sun's radius
pars.Astro.T_sun        = 5778;            % [k]           Temperature of the Sun

% SRP parameters
pars.SRP.eta   = 1; % SRP efficiency
pars.SRP.d0    = [-2.54699305407e7;1.329307959473e8;5.762458661619e7]; % position vector from sun to primary in km
pars.SRP.const = pars.Astro.R_sun^2*pars.Astro.boltzman*pars.Astro.T_sun^4/pars.Astro.speedOfLight; % constant for solar radiation pressure in W*s/km = kg*km/s^2

% third body perturbations parameters
pars.TBP.GM = 1.3271244004127e11; % GMs of third bodies in km^3/s^2
pars.TBP.d0 = [-2.54699305407e7;1.329307959473e8;5.762458661619e7]; % position vectors from third bodies to primary in km
 
% s/c states
xsc0 = [
    [757.7;
     5222.607;
     4851.5;
     2.213;
     4.67834;
     -5.3713],...
    [-6120.7347;
    20606.4449;
    -125.8035;
    -4.130255006;
    -1.22457538;
    -0.039121849]
    ]; % s/c initial states in km and km/s

% output measurement files
file1 = 'Project1_SC%d_Measurements_Raw_Nominal.txt';
file2 = 'Project1_SC%d_Measurements_Noisy_Nominal.txt';
file3 = 'Project1_SC%d_Trajectory_Data_Nominal_ObsTimes.txt';
file4 = 'Project1_SC%d_Trajectory_Data_Nominal_AllTimes.txt';

%% For quick access
t0         = pars.Body.t0;
Wt0        = pars.Body.Wt0;
Wdot       = pars.Body.Wdot;
RA         = pars.Body.RA;
DEC        = pars.Body.DEC;
x_stations = pars.Body.Stations.State;
IDs        = pars.Body.Stations.IDs;
h_hat      = pars.Body.h_hat; % angular velocity vector of primary in rad/s 
el_cutoff  = 10*pi/180; % Elevation cutoff below which we don't create measurements
