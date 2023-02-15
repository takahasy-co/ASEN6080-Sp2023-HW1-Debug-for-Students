function out = tilde(vec)

r1  = vec(1);
r2  = vec(2);
r3  = vec(3);
out = [0,-r3,r2;
       r3,0,-r1;
       -r2,r1,0];
