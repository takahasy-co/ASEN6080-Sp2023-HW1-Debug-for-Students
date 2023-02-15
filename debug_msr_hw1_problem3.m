% This script validates the analytical partials by finite differencing
% of the observables

addpath('lib')
addpath('lib/util')
addpath('lib/propagator')
addpath('lib/measurements')
addpath('lib/inputs')

clear
clc
format long g

loadDebugHW1;

r_N = xsc0(1:3,1);
v_N = xsc0(4:6,1);
x_N = [r_N;v_N];

t0         = pars.Body.t0;
Wt0        = pars.Body.Wt0;
Wdot       = pars.Body.Wdot;
RA         = pars.Body.RA;
DEC        = pars.Body.DEC;
x_stations = pars.Body.Stations.State;
IDs        = pars.Body.Stations.IDs;
h_hat      = pars.Body.h_hat; % angular velocity vector of primary in rad/s

tf = 7*24*3600;
dt = 60;
t  = t0+(0:dt:tf);
ii = 100;

% Compute the rotation matrix between the inertial and body frame
Wt  = wrapTo2Pi(Wt0 + Wdot*(t(ii)-t0)); % position of prime meridian in rad
BN  = BodyFrame(Wt,DEC,RA);             % inertial to body frame rotation matrix
NB  = BN';                              % body frame to inertial frame rotation matrix

stnIdx      = 1;
x_station_B = x_stations(:,stnIdx);
r_station_B = x_station_B(1:3);
v_station_B = x_station_B(4:6); % Should be [0;0;0]
r_station_N = NB*r_station_B;
v_station_N = NB*( v_station_B + cross(h_hat,r_station_B) );
x_station_N = [r_station_N; v_station_N];

[rho,rhodot] = Range_RangeRate(x_N,x_station_N);

rho_vec      = r_N - r_station_N;
rho_hat      = rho_vec/rho;

d_rho_d_r_sc = rho_hat';
d_rho_d_v_sc = zeros(1,3);
d_rho_d_x_sc = [d_rho_d_r_sc, d_rho_d_v_sc];

d_rhodot_d_r_sc = 1/rho*(v_N - v_station_N)'*(eye(3) - rho_hat*rho_hat');
d_rhodot_d_v_sc = rho_hat';
d_rhodot_d_x_sc = [d_rhodot_d_r_sc, d_rhodot_d_v_sc];

d_rho_d_Rs     = -rho_hat'*NB;
d_rho_dot_d_Rs = -1/rho * ( ( v_N - v_station_N )' * ( eye(3) - rho_hat*rho_hat' ) * NB + rho_vec'*(NB)*tilde(h_hat) );

delta      = 1e-8;
par_sc     = zeros(1,6);
pardot_sc  = zeros(1,6);
par_stn    = zeros(1,3);
pardot_stn = zeros(1,3);

% rho and rhodot partials w.r.t. s/c state

for ii = 1:6

    delta_x_sc     = zeros(6,1);
    delta_x_sc(ii) = delta;
    x_N_p          = x_N + delta_x_sc;
    
    [rho_p,rhodot_p] = Range_RangeRate(x_N_p,x_station_N);
    
    par_sc(ii)    = (rho_p - rho)/delta;
    pardot_sc(ii) = (rhodot_p - rhodot)/delta;
    
end

par_ratio_rho_sc    = d_rho_d_x_sc./par_sc;
par_ratio_rhodot_sc = d_rhodot_d_x_sc./pardot_sc;

d_rho_d_x_sc
par_sc
par_ratio_rho_sc

d_rhodot_d_x_sc
pardot_sc
par_ratio_rhodot_sc

% rho and rhodot partials w.r.t. body frame station location

for ii = 1:3

    x_station_B   = x_stations(:,stnIdx);
    delta_Rs      = zeros(3,1);
    delta_Rs(ii)  = delta;
    r_station_B_p = x_station_B(1:3) + delta_Rs;
    v_station_B   = x_station_B(4:6); % Should be [0;0;0]
    r_station_N_p = NB*r_station_B_p;
    v_station_N_p = NB*( v_station_B + cross(h_hat,r_station_B_p) );
    x_station_N_p = [r_station_N_p; v_station_N_p];
    
    [rho_p,rhodot_p] = Range_RangeRate(x_N,x_station_N_p);
    
    par_stn(ii)    = (rho_p - rho)/delta;
    pardot_stn(ii) = (rhodot_p - rhodot)/delta;
    
end

par_ratio_rho_stn    = d_rho_d_Rs./par_stn;
par_ratio_rhodot_stn = d_rho_dot_d_Rs./pardot_stn;

d_rho_d_Rs
par_stn
par_ratio_rho_stn

d_rho_dot_d_Rs
pardot_stn
par_ratio_rhodot_stn
