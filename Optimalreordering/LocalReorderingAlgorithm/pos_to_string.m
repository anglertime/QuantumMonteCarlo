function s = pos_to_string(x)

global dim

if dim == 1
    s = sprintf('[ %f ]',x(1));
elseif dim == 2
    s = sprintf('[ %f %f]',x(1),x(2));
elseif dim == 3
    s = sprintf('[ %f %f %f ]',x(1),x(2),x(3));
end