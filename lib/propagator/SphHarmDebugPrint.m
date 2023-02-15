% This script will integrate the equations of motion for state
% under the influence of a spherical harmonic expansion of a gravity field,
% written for ASEN 6080 Spring 2023

function [Acce_body,Acce_inertial,A_Acce_Pos_B,A_Acce_Pos_N,A_Acce_C_B,A_Acce_C_N,A_Acce_S_B,A_Acce_S_N] = SphHarmDebugPrint(t,x,pars)

% GM
Rref   = pars.Rref; % reference radius
GM     = pars.GM;   % GM

% SH
n_degree = pars.SH.degree;
C        = pars.SH.C;
S        = pars.SH.S;

% Pole and body frame
RA   = pars.RA;  % right ascension b/t inertial and body frame in rad
DEC  = pars.DEC; % declination b/t inertial and body frame in rad
Wt0  = pars.Wt0; % body's rotation rate rad/s
Wdot = pars.Wdot;
t0   = pars.t0;
Wt   = wrapTo2Pi(Wt0 + Wdot*(t-t0));
BN   = BodyFrame(Wt,DEC,RA); % inertial to body frame
NB   = BN';

% convert inertial pos to body frame pos
r_N = x(1:3); % position in inertial frame
r_B = BN*r_N; % position in body frame

% Velocity
v_N = x(4:6);

% Acceleration
[~, Acce_body, A_Acce_Pos_B, A_Acce_C_B, A_Acce_S_B] = ComputePotentialAcceSTMConsider_mex(n_degree, n_degree, n_degree, Rref, GM, r_B, C, S);
Acce_inertial = NB*Acce_body;
A_Acce_Pos_N  = NB*A_Acce_Pos_B*BN;
A_Acce_C_N    = NB*A_Acce_C_B;
A_Acce_S_N    = NB*A_Acce_S_B;

end