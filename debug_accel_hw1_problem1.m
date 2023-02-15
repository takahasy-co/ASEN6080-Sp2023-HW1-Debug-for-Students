% This script checks for the acceleration computation consistency

addpath('lib')
addpath('lib/util')
addpath('lib/propagator')
addpath('lib/measurements')
addpath('inputs')

clear
clc
format long g

loadDebugHW1;

% Time

t = 100; % seconds past epoch

% GM
Rref   = pars.Body.Rref; % reference radius
GM     = pars.Body.GM;   % GM

% SH
n_degree = pars.Body.SH.degree;

% Pole and body frame
RA   = pars.Body.RA; % right ascension b/t inertial and body frame in rad
DEC  = pars.Body.DEC; % declination b/t inertial and body frame in rad
Wt0  = pars.Body.Wt0; % body's rotation rate rad/s
Wdot = pars.Body.Wdot;
t0   = pars.Body.t0;
Wt   = wrapTo2Pi(Wt0 + Wdot*(t-t0));
BN   = BodyFrame(Wt,DEC,RA); % inertial to body frame
NB   = BN';

% convert inertial pos to body frame pos
r_N = xsc0(1:3); % position in inertial frame
r_B = BN*r_N; % position in body frame

% Acceleration
[U, Acce_body, A_Acce_Pos, A_Acce_C, A_Acce_S] = ComputePotentialAcceSTMConsider_mex(n_degree, n_degree, n_degree, Rref, GM, r_B, pars.Body.SH.C, pars.Body.SH.S);

Acce_inertial       = NB*Acce_body;
A_Acce_Pos_inertial = NB*A_Acce_Pos*BN;

A_Acce_C_inertial = NB*A_Acce_C;
A_Acce_S_inertial = NB*A_Acce_S;

acce_from_CS_body   = Acce_body - (-GM/norm(r_B)^3*r_B);        % Total acceleration minus two-body
acce_from_CS_body_2 = [A_Acce_C, A_Acce_S]*pars.Body.SH.CS_vec; % Same as above but computed in the vector format

U
Acce_body
Acce_inertial

A_Acce_Pos
A_Acce_Pos_inertial

A_Acce_C
A_Acce_C_inertial

A_Acce_S
A_Acce_S_inertial

acce_from_CS_body
acce_from_CS_body_2
