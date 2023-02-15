% This script will integrate the equations of motion for state
% under the influence of a spherical harmonic expansion of a gravity field,
% written for ASEN 6080 Spring 2023

function dx = SphHarmWithPhi(t,x,pars,prop)

% GM
Rref   = pars.Rref; % reference radius
GM     = pars.GM;   % GM

% Update GM if it's in the state vector
for ii = 1:length(prop.biasNames)

    biasIdx = prop.biasIdx(ii);
    
    if strcmp(prop.biasNames{ii},'GM')
        
        GM = x(biasIdx);
        
    end % For if
    
end % For ii

% SH
n_degree = pars.SH.degree;
C        = pars.SH.C;
S        = pars.SH.S;

% Update C and S harmonics matrices if it's in the state vector
for ii = 1:length(prop.biasNames)

    biasIdx = prop.biasIdx(ii);
    
    if strcmp(prop.biasNames{ii}(1),'C')
        
        nn             = str2double(prop.biasNames{ii}(2)); % This only works for < 10 degree harmonics
        mm             = str2double(prop.biasNames{ii}(3)); % This only works for < 10 degree harmonics
        n_idx          = nn + 1;
        m_idx          = mm + 1;
        c_now          = x(biasIdx);
        C(n_idx,m_idx) = c_now;
        
    elseif strcmp(prop.biasNames{ii}(1),'S') % This only works for < 10 degree harmonics
        
        nn             = str2double(prop.biasNames{ii}(2)); % This only works for < 10 degree harmonics
        mm             = str2double(prop.biasNames{ii}(3)); % This only works for < 10 degree harmonics
        n_idx          = nn + 1;
        m_idx          = mm + 1;
        s_now          = x(biasIdx);
        S(n_idx,m_idx) = s_now;
        
    end % For if
    
end % For ii

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
[~, Acce_body, A_Acce_Pos, A_Acce_C, A_Acce_S] = ComputePotentialAcceSTMConsider_mex(n_degree, n_degree, n_degree, Rref, GM, r_B, C, S);
Acce_inertial = NB*Acce_body;

% STM computation

numParams = prop.numParams;
Phi_dot   = reshape(zeros(numParams,numParams),numParams^2,1);

if prop.partials
    
    A          = zeros(numParams,numParams);
    A(1:3,4:6) = eye(3);
    A(4:6,1:3) = NB*A_Acce_Pos*BN;
    
    for ii = 1:length(prop.biasNames)
        
        biasIdx = prop.biasIdx(ii);
        
        if strcmp(prop.biasNames{ii},'GM')
            
            A(4:6,biasIdx) = Acce_inertial/GM;
            
        elseif strcmp(prop.biasNames{ii}(1),'C')
            
            nn = str2double(prop.biasNames{ii}(2)); % This only works for < 10 degree harmonics
            mm = str2double(prop.biasNames{ii}(3)); % This only works for < 10 degree harmonics
            
            % printf('%i %i %s\n',nn,mm,prop.biasNames{ii})
            
            c_idx          = nn*(nn+1)/2+mm+1-1; % -1 for C00
            A(4:6,biasIdx) = NB*A_Acce_C(:,c_idx);
            
        elseif strcmp(prop.biasNames{ii}(1),'S') % This only works for < 10 degree harmonics
            
            nn = str2double(prop.biasNames{ii}(2)); % This only works for < 10 degree harmonics
            mm = str2double(prop.biasNames{ii}(3)); % This only works for < 10 degree harmonics
            
            s_idx = nn*(nn-1)/2+mm;
            
            A(4:6,biasIdx) = NB*A_Acce_S(:,s_idx);
            
        end % For if
        
    end % For ii
    
    Phi       = reshape(x(numParams+1:end),numParams,numParams);
    Phi_dot   = reshape(A*Phi,numParams^2,1);
    
end % For if

DdynamicDt = zeros(prop.numDynamic,1);
DbiasDt    = zeros(prop.numBias,1);
DstochDt   = zeros(prop.numStoch,1);

% Time derivative
dx  = [v_N;Acce_inertial;DdynamicDt;DbiasDt;DstochDt;Phi_dot];

end