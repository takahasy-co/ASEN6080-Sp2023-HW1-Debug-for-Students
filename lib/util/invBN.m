function [gamma,beta,alpha] = invBN(BN)

alpha = atan2(BN(3,1),-BN(3,2));
beta  = acos(BN(3,3));
gamma = atan2(BN(1,3),BN(2,3));