function x = put_in_box (y)

% Shifts any particles from outside orbital cube to inside cube by
% utilizing the periodicity

global Box

% x = y-floor((y + 1/4*sideLength)/Box(1))*Box(1);
%2013-04-23 Ming Li
%Change the cube of particles into [0, L)^3
x = y-floor(y/Box(1))*Box(1);