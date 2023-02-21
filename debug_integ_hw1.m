% This script outputs all sorts of dynamical/measurement parameters
% at various times for students' debugging.

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
loadPropDebugIntegHW1

% create trajectory
[t,y] = ode45(@(t,x) sc_EOM_Phi(t,x,pars,prop),T,x0,prop.opts);

numTime                = length(t);

trajectory_data        = zeros(numTime,7);
acce_B_data            = zeros(numTime,4);
acce_N_data            = zeros(numTime,4);
A_Acce_Pos_B_data      = zeros(numTime,10);
A_Acce_Pos_N_data      = zeros(numTime,10);
A_Acce_C_B_data        = zeros(numTime,3*(pars.Body.SH.num_C-1)+1); % -1 for C00, +1 for time
A_Acce_C_N_data        = zeros(numTime,3*(pars.Body.SH.num_C-1)+1); % -1 for C00, +1 for time
A_Acce_S_B_data        = zeros(numTime,3*(pars.Body.SH.num_S)+1);
A_Acce_S_N_data        = zeros(numTime,3*(pars.Body.SH.num_S)+1);

Phi_data               = zeros(numTime,prop.numParams^2+1);

bn_data                = zeros(numTime,10);
measurement_data_raw   = zeros(3*numTime,4);
measurement_data_raw_2 = zeros(3*numTime,4);
azel_data              = zeros(3*numTime,4);

% loop over time
for ii = 1:length(t)
        
    % Compute the rotation matrix between the inertial and body frame
    Wt  = wrapTo2Pi(Wt0 + Wdot*(t(ii)-t0)); % position of prime meridian in rad
    BN  = BodyFrame(Wt,DEC,RA);             % inertial to body frame rotation matrix
    NB  = BN';                              % body frame to inertial frame rotation matrix
    
    % BN data
    bn_data(ii,:) = [t(ii),reshape(BN,1,9)];
    
    % trajectory information
    trajectory_data(ii,:) = [t(ii),y(ii,1:6)];
    Phi_data(ii,:)        = [t(ii),y(ii,prop.numParams+1:end)];
    
    % Acceleration data
    [Acce_B,Acce_N,A_Acce_Pos_B,A_Acce_Pos_N,A_Acce_C_B,A_Acce_C_N,A_Acce_S_B,A_Acce_S_N] = SphHarmDebugPrint(t(ii),y(ii,1:6)',pars.Body);
    acce_B_data(ii,:)       = [t(ii),reshape(Acce_B,1,3)];
    acce_N_data(ii,:)       = [t(ii),reshape(Acce_N,1,3)];
    A_Acce_Pos_B_data(ii,:) = [t(ii),reshape(A_Acce_Pos_B,1,9)];
    A_Acce_Pos_N_data(ii,:) = [t(ii),reshape(A_Acce_Pos_N,1,9)];
    A_Acce_C_B_data(ii,:)   = [t(ii),reshape(A_Acce_C_B,1,3*(pars.Body.SH.num_C-1))];
    A_Acce_C_N_data(ii,:)   = [t(ii),reshape(A_Acce_C_N,1,3*(pars.Body.SH.num_C-1))];
    A_Acce_S_B_data(ii,:)   = [t(ii),reshape(A_Acce_S_B,1,3*(pars.Body.SH.num_S))];
    A_Acce_S_N_data(ii,:)   = [t(ii),reshape(A_Acce_S_N,1,3*(pars.Body.SH.num_S))];
    
    for stnIdx = 1:length(pars.Body.Stations.IDs)
        
        % Station locations in inertial and body frames
        x_station_B = x_stations(:,stnIdx);
        r_station_B = x_station_B(1:3);
        v_station_B = x_station_B(4:6); % Should be [0;0;0]
        r_station_N = NB*r_station_B;
        v_station_N = NB*( v_station_B + cross(h_hat,r_station_B) );
        x_station_N = [r_station_N; v_station_N];
        
        % s/c state in the inertial frame
        r_sc_N = y(ii,1:3)';    % s/c position in inertial frame
        v_sc_N = y(ii,4:6)';    % s/c velocity in inertial frame
        x_sc_N = [r_sc_N; v_sc_N];
        
        [az,el]                      = azel(r_sc_N,BN,r_station_B);
        azel_data(3*(ii-1)+stnIdx,:) = [t(ii),IDs(stnIdx), az*180/pi, el*180/pi];

        % compute measurement
        [rho,rhodot]                            = Range_RangeRate(x_sc_N,x_station_N);
        measurement_data_raw(3*(ii-1)+stnIdx,:) = [t(ii),IDs(stnIdx),rho,rhodot];

        % This is what you get when you feed the inertial quantities expressed in body frame
        % Should be identical to the inertial formulation
        x_sc_N_2                                  = zeros(6,1);
        x_station_N_2                             = zeros(6,1);
        x_sc_N_2(1:3)                             = BN*r_sc_N;
        x_sc_N_2(4:6)                             = BN*v_sc_N;
        x_station_N_2(1:3)                        = r_station_B;
        x_station_N_2(4:6)                        = cross(h_hat,r_station_B);
        [rho_2,rhodot_2]                          = Range_RangeRate(x_sc_N_2,x_station_N_2);
        measurement_data_raw_2(3*(ii-1)+stnIdx,:) = [t(ii),IDs(stnIdx),rho_2,rhodot_2];
        
    end % For stnIdx

end % For ii

% save data
formatSpec = '%1.6e %3d %+1.15e %+1.15e\n';
Write_Data('hw1_debug_outputs/raw_obs_no_noise_debug.txt',measurement_data_raw,formatSpec)
Write_Data('hw1_debug_outputs/raw_obs_no_noise_debug_2.txt',measurement_data_raw_2,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,6),'\n');
Write_Data('hw1_debug_outputs/traj_data_debug.txt',trajectory_data,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,prop.numParams^2),'\n');
Write_Data('hw1_debug_outputs/phi_debug.txt',Phi_data,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,9),'\n');
Write_Data('hw1_debug_outputs/bn_matrix_debug.txt',bn_data,formatSpec)

formatSpec = '%1.6e %3d %+1.15e %+1.15e\n';
Write_Data('hw1_debug_outputs/az_el_debug.txt',azel_data,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,3),'\n');
Write_Data('hw1_debug_outputs/accel_B_debug.txt',acce_B_data,formatSpec)
Write_Data('hw1_debug_outputs/accel_N_debug.txt',acce_N_data,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,9),'\n');
Write_Data('hw1_debug_outputs/A_matrix_dAcce_dPos_B_debug.txt',A_Acce_Pos_B_data,formatSpec)
Write_Data('hw1_debug_outputs/A_matrix_dAcce_dPos_N_debug.txt',A_Acce_Pos_N_data,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,3*(pars.Body.SH.num_C-1)),'\n');
Write_Data('hw1_debug_outputs/A_matrix_dAcce_dC_B_debug.txt',A_Acce_C_B_data,formatSpec)
Write_Data('hw1_debug_outputs/A_matrix_dAcce_dC_N_debug.txt',A_Acce_C_N_data,formatSpec)

formatSpec = strcat('%1.6e',repmat(' %+1.15e',1,3*pars.Body.SH.num_S),'\n');
Write_Data('hw1_debug_outputs/A_matrix_dAcce_dS_B_debug.txt',A_Acce_S_B_data,formatSpec)
Write_Data('hw1_debug_outputs/A_matrix_dAcce_dS_N_debug.txt',A_Acce_S_N_data,formatSpec)

plot3(y(:,1),y(:,2),y(:,3)); axis equal

%% EOM
function dx = sc_EOM_Phi(t,x,pars,prop)

    % GM + SH
    dx = accel_hw1(t,x,pars.Body,prop);
    
end % For function
