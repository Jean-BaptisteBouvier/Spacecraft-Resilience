%%% Function setting up all the parameters for the main code


function params = parameters_setup(dt, failure, mass, thrust, transfer_time, waypoints)

%%% Waypoint following Chapter 7 of Ortolano's thesis
params = CW_RPO_TestCondition(125102990);
params.Omega = 0.001060922896439; % [rad/s]
params.simTimeHours = transfer_time; % [hours]
params.vInit = [0; 0];
params.vFinal = [0; 0];
params.mass = mass; % [kg]
params.thrust_factor = thrust/params.mass; % = 1.5e-4 [m/s^2] thrust-to-weight ratio
params.dt = dt;
params.waypoints = waypoints;
params.V_exit = 1660*9.81; % [m/s] = Isp [s] x g [m/s^2] exit velocity of ions in PPS-1350
params.failure = failure; % id of the malfunctioning thruster in {1,2,3,4,5}. Only resilient to no. 4

%%% Matrices of the system
B_bar = params.thrust_factor*[zeros(2,5); 1, 1, -1, -sqrt(2), -1; 1, -1, -1, 0, 1];

B = B_bar; B(:, failure) = [];
params.matrix_B = B;
C = B_bar(:, failure);
params.matrix_C = C;
A = zeros(4,4);
A(3,1) = 3*params.Omega^2; A(1,3) = 1; A(2,4) = 1;
A(3,4) = 2*params.Omega; A(4,3) = -2*params.Omega;
params.matrix_A = A;


end