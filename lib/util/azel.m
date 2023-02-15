function [az,el] = azel(R_N,BN,Rs_B)

B_ENU = ENU(Rs_B);
ENU_B = B_ENU';

NB    = BN';
Rs_N  = NB*Rs_B;

dR_ENU = ENU_B*BN*(R_N-Rs_N);

el = asin(dR_ENU(3)/norm(dR_ENU));
az = atan2(dR_ENU(2),dR_ENU(1));