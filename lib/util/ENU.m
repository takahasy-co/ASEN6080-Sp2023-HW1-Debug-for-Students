function B_ENU = ENU(Rs_B)

Up    = Rs_B./norm(Rs_B);
if all(Up == [0;0;1])
   East  = cross([0;-1;0],Up);
else
   East  = cross([0;0;1],Up);
end
East  = East./norm(East);
North = cross(Up,East);
North = North./norm(North);

B_ENU   = [East North Up];
