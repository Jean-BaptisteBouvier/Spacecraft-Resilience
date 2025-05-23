%%% Main file for the resilient trajectory tracking to observe a target satellite
%%%
%%% Authors: Jean-Baptiste Bouvier and Himmat Panag.


clearvars
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%% User inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delay = 3; % [s]
dt = 0.2; % [s] time step must be smaller than the delay

is_Lipschitz_w = true; % [true]: w is Lipschitz   [false]: w is bang-bang
w_max = 1; % maximal value of w, must be <= 1
L = 0.1; % Lipschitz constant of w
num_bang_per_hour = 10; % average number of bangs per hour when w is bang-bang

mass = 600; % [kg] spacecraft mass
thrust = 90e-3; % [N] max thrust of each thruster (90 mN for the PPS-1350)
failure = 4; % id of the malfunctioning thruster in {1,2,3,4,5}.  Only resilient to no. 4
transfer_time = 1.5; % [hours] transfer time between each waypoints
waypoints = [0,200; 0,80; 80,0; 0,-80; -80,0; 0,80].*1e-3; % [km] waypoints positions




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loss of control authority over thruster no. ' + string(failure))
disp('Actuation delay = '+ string(delay) + 's       dt = ' +string(dt)+'s    Transfer time = ' + string(60*transfer_time) + 'min')

%%% Add path to the folders of functions and data
addpath 'functions'
addpath 'data'
%%% Setting up all the parameters
params = parameters_setup(dt, failure, mass, thrust, transfer_time, waypoints);
N = transfer_time*60*60/dt; % number of time steps per transfer
nb_transfers = length(waypoints(:,1))-1; % number of transfers
p = length(params.matrix_C(1,:)); % number of malfunctioning thrusters
m = length(params.matrix_B(1,:)); % number of controlled thrusters


%%% Generating the reference trajectory if not already done
filename = 'data/ref_traj_' + string(60*transfer_time) + 'min_dt=' + string(dt) + '.mat';
if isfile(filename)
    load(filename)
else
    [X_ref, U_ref] = reference_trajectory(params);
end
%%% Reference trajectory analysis
reference_trajectory_analysis(params, X_ref, U_ref)

%%% Generating the undesirable inputs if not already done
W = undesirable_input(is_Lipschitz_w, w_max, L, num_bang_per_hour, params);

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

