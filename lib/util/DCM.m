%
% Yu Takahashi modified the original script by David Lujan.
% The order of matrix operation is from right to left.
%
% DCM.m creates a DCM to produce a sequence of rotations using Euler angles
%   to transform from one frame to another.
%
% function y = DCM(type,angles)
%
% Input:
%   type: n-element array with the order of rotations to perform.
%           ex: type = [1,2,3,2,3] is a 3-2-3-2-1 rotation
%   angles: n-element array with the angles in radians to use for each axis
%       rotation.
%
% Output:
%   y: [3x3] proper orthogonal matrix
%
% Dependencies: none
%
% Author: David Lujan, david.lujan@colorado.edu

function y = DCM(type,angles)
n = length(type);
y = eye(3);
for k = n:-1:1
    a = angles(k);
    %fprintf('%d-th input, %d-th rotation, axis = %d, angle = %f\n',k,n-k+1,type(k),a*180/pi)
    if type(k) == 1
        M = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
    elseif type(k) == 2
        M = [cos(a) 0 -sin(a);0 1 0;sin(a) 0 cos(a)];
    elseif type(k) == 3
        M = [cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1];
    else
        error('Incorrect rotation type. Must be a 1, 2, or 3')
    end
    y = M*y;
end
end