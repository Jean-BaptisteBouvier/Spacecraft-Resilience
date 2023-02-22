%%% Main file for the resilient trajectory tracking to observe a target satellite
%%%
%%% Authors: Jean-Baptiste Bouvier and Himmat Panag.


clearvars
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%% User inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delay = 1; % [s]
dt = 1; % [s] time step must be smaller than the delay

is_Lipschitz_w = true; % [true]: w is Lipschitz   [false]: w is bang-bang
w_max = 1e-2; % maximal value of w, must be < 1
dwdt_max = w_max*1e-3; % max time derivative of w
num_bang_per_hour = 10; % average number of bangs per hour when w is bang-bang

mass = 600; % [kg] spacecraft mass
thrust = 90e-3; % [N] max thrust of each thruster (90 mN for the PPS-1350)
failure = 4; % id of the malfunctioning thruster in {1,2,3,4,5}.  Only resilient to no. 4
transfer_time = 1.5; % [hours] transfer time between each waypoints
waypoints = [0,200; 0,80; 80,0; 0,-80; -80,0; 0,80].*1e-3; % [km] waypoints positions




%%%%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loss of control authority over thruster no. ' + string(failure))
disp('Actuation delay = '+ string(delay) + 's     Transfer time = ' + string(60*transfer_time) + 'min')

%%% Add path to the folder of functions
addpath 'functions'
%%% Setting up all the parameters
params = parameters_setup(dt, failure, mass, thrust, transfer_time, waypoints);
N = transfer_time*60*60/dt;
p = length(params.matrix_C(1,:)); % number of malfunctioning thrusters
m = length(params.matrix_B(1,:)); % number of controlled thrusters
nb_transfers = length(waypoints(:,1))-1;


%%% Generating the reference trajectory if not already done
filename = 'ref_traj_' + string(60*transfer_time) + 'min.mat';
if isfile(filename)
    load(filename)
else
    [X_ref, U_ref] = reference_trajectory(params);
end
%%% Reference trajectory analysis
reference_trajectory_analysis(params, X_ref, U_ref)


%%% Generating the undesirable inputs if not already done
W = zeros(p, N*nb_transfers);
if is_Lipschitz_w 
	W_Lip_filename = 'W_Lip_' + string(60*transfer_time) + 'min.mat';
	disp('w is Lipschitz and w_max = ' + string(w_max) + 'm/s^2')
	if isfile(W_Lip_filename)
		load(W_Lip_filename);
    else
		W(:,1) = w_max*rand(p,1);
		for i = 2:N*nb_transfers
			W(:,i) = W(:,i-1) + params.dt*dwdt_max*2*(rand(p,1)-0.5);
			W(:,i) = ((W(:,i) <= w_max).*(W(:,i) >= 0)).*W(:,i) + (W(:,i) > w_max)*w_max;
		end
		save(W_Lip_filename, 'W');
	end
else % w is bang-bang
	W_bang_filename = 'W_bang_' + string(60*transfer_time) + 'min.mat';
	disp('w is bang-bang and w_max = ' + string(w_max) + 'm/s^2')
	if isfile(W_bang_filename)
		load(W_bang_filename);
    else
        bang_probability = num_bang_per_hour * dt/3600;
		W(:,1) = w_max*ones(p,1);
		for i = 2:N*nb_transfers
			if rand > 1-bang_probability
				W(:,i) = (W(:,i-1) == zeros(p,1))*w_max;
			else
				W(:,i) = W(:,i-1);
			end
		end
		save(W_bang_filename, 'W');
	end
end


%%% Generating the tracking trajectory for each transfer between waypoints
X_Lechappe = zeros(4,N*nb_transfers); U_Lechappe = zeros(4,N*nb_transfers);

for transferNum = 1:nb_transfers
    interval_id = 1+N*(transferNum-1):N*transferNum;

    %%% Selecting the transfer start point
    if transferNum == 1
        x0 = X_ref(:,1);
    else
        x0 = x_L(:,end);
    end    

    %%% Trajectory tracking
    [x_L, u_L] = Lechappe_tracking(x0, X_ref(:, interval_id), U_ref(:, interval_id), W(:, interval_id), delay, params);
    
    %%% Storing data
    X_Lechappe(:,interval_id) = x_L;
    U_Lechappe(:,interval_id) = u_L;
end

%%% Tracking trajectory analysis
trajectory_analysis(params, X_ref, U_ref, W, X_Lechappe, U_Lechappe)

