% Body parameters
pars.Body.GM        = 4.19300435e5;           % Gravitational parameter km^3/s^2
pars.Body.Rref      = 6378.136;               % reference radius for spherical harmonics in km
pars.Body.t0        = 0;                      % epoch for prime meridian
pars.Body.Wt0       = 10.5*pi/180;            % location of prime meridian at epoch in rad
pars.Body.Wdot      = 2*pi/(2023);            % rotation rate of prime meridian in rad/s
pars.Body.h_hat     = pars.Body.Wdot*[0;0;1]; % angular velocity vector of body fixed frame relative to the inertial frame in rad/s 
pars.Body.RA        = 52.9*pi/180;            % right ascension of Body pole in rad
pars.Body.DEC       = 63.5*pi/180;            % declination of Body pole in rad
                    
pars.Body.SH.degree = 3; % The degree and order of the harmonics.

pars.Body.SH.C      = zeros(pars.Body.SH.degree+1,pars.Body.SH.degree+1);
pars.Body.SH.C(1,1) = 1.0;
pars.Body.SH.C(3,1) = -2.0828335e-3; % C20
pars.Body.SH.C(3,2) = -8.5898e-5;    % C21
pars.Body.SH.C(3,3) = 5.57457e-3;    % C22
pars.Body.SH.C(4,1) = 4.518373e-3;   % C30
pars.Body.SH.C(4,2) = -3.516418e-4;  % C31
pars.Body.SH.C(4,3) = -8.87191e-3;   % C32
pars.Body.SH.C(4,4) = 7.15129e-5;    % C33

pars.Body.SH.S      = zeros(pars.Body.SH.degree+1,pars.Body.SH.degree+1);
pars.Body.SH.S(3,2) = 1.80326e-4;    % S21
pars.Body.SH.S(3,3) = -9.03867e-7;   % S22
pars.Body.SH.S(4,2) = 8.1412e-5;     % S31
pars.Body.SH.S(4,3) = -5.111891e-5;  % S32
pars.Body.SH.S(4,4) = 4.82403e-4;    % S33

pars.Body.SH.num_C  = (pars.Body.SH.degree+1)*(pars.Body.SH.degree+2)/2;
pars.Body.SH.num_S  = pars.Body.SH.degree*(pars.Body.SH.degree+1)/2;
pars.Body.SH.num_CS = pars.Body.SH.num_C + pars.Body.SH.num_S;
pars.Body.SH.C_vec  = zeros(pars.Body.SH.num_C-1,1); % -1 for C00
pars.Body.SH.S_vec  = zeros(pars.Body.SH.num_S,1);

for nn = 1:pars.Body.SH.degree
    
    n_idx = nn + 1;
    
    for mm = 0:nn
    
        m_idx = mm + 1;
        Cnm   = pars.Body.SH.C(n_idx,m_idx);
        
        c_idx = nn*(nn+1)/2+mm+1-1; % -1 for C00
        
        pars.Body.SH.C_vec(c_idx) = pars.Body.SH.C(n_idx,m_idx);
        
        if mm > 0
        
            s_idx = nn*(nn-1)/2+mm;
            pars.Body.SH.S_vec(s_idx) = pars.Body.SH.S(n_idx,m_idx);
        
        end % For if
        
    end % For mm
            
end % For ii

pars.Body.SH.CS_vec = [pars.Body.SH.C_vec; pars.Body.SH.S_vec];

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

% s/c states
xsc0 = [
    757.7;
     5222.607;
     4851.5;
     2.213;
     4.67834;
     -5.3713
    ]; % s/c initial states in km and km/s

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
