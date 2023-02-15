function BN = BodyFrame(Wt,DEC,RA)

if DEC == pi/2
   BN = DCM([3,1,3],[Wt,(pi/2-DEC),RA]);
else
   BN = DCM([3,1,3],[Wt,(pi/2-DEC),(pi/2+RA)]);
end % For if