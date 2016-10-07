function x = convert_dist( y )
% when the distance of particle and orbital is larger than one half of the
% Box length, convert it.

% This is only for one axis.

global Box dim datafile

% L = period;
% x = y-floor(y/L+0.5)*L; %% Map R into [-0.5L, 0.5L).

% 2012-Nov-5 Change the difference into positive distance. Only works if
% % distance is small than L.
dif = abs(y);
x = y;
for i = 1:dim
    if (dif(i) > 0.5*Box(i))
        x(i) = Box(i) - dif(i);
    else
        x(i) = dif(i);
    end
    if (dif(i) >= Box(i))
        fprintf(datafile, 'Error! The distance is larger than or equal to L \n');
    end
end
end