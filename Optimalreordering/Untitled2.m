global constK N_length sideln max_dist period

N_length = 40; dim = 1;

constK = 0.5;
[Pconfig_ideal, Oconfig, period, ~, sideln, max_dist] = ...
    read_configuration_max_dis (1024, dim, N_length);

% perturb = 0.3*sideln;
% temp = Pconfig_ideal + 0.3*sideln;
Pconfig_ideal = put_in_box(Pconfig_ideal + 0.3*sideln);

Pconfig_ideal_moved = Pconfig_ideal;
Pconfig_ideal_moved(35) = Pconfig_ideal_moved(35) - 15.1*sideln;
Pconfig_ideal_moved(34) = Pconfig_ideal_moved(34) - 16.4*sideln;
Pconfig_ideal_moved(30) = Pconfig_ideal_moved(30) - 12.1*sideln;
Pconfig_ideal_moved = put_in_box(Pconfig_ideal_moved);
% Pconfig_ideal_moved(26) = Pconfig_ideal_moved(26) - 14.7*sideln;

A_ideal = compute_matrix (Pconfig_ideal, Oconfig);
A_shifted = compute_matrix (Pconfig_ideal_moved, Oconfig);