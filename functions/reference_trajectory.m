%%% Generating the reference trajectory

%%% Based on the method from the paper:
%%% "Autonomous Optimal Trajectory Planning for Orbital Rendezvous,
%%% Satellite Inspection, and Final Approach Based on Convex Optimization"
%%% from Nicholas Ortolano, David K. Geller, and Aaron Avery.
%%% Non rotating target. Optimal trajectories for the Hill-Clohessy-Wiltshire 
%%% model in cartesian coordinates using SOCPs. 

function [X_ref, U_ref] = reference_trajectory(params)

dt = params.dt; % desired time step
N = params.transfer_time*60*60/dt; % number of steps per transfer
nb_transfers = length(params.waypoints(:,1))-1;
[n, m] = size(params.matrix_B);
X_ref = zeros(n, N*nb_transfers);
U_ref = zeros(m, N*nb_transfers);


if dt > 1 % for large dt calculations can be done in one iteration
    params.numSteps = N;
    for transferNum = 1:nb_transfers
        %%% Define start and end point of the transfer
        params.rInit = params.waypoints(transferNum,:)'; % x is r-bar, y is v-bar, z is cross track        
        params.rFinal = params.waypoints(transferNum+1,:)';
        %%% Calculate trajectory without the spacecraft rotation
        [~, x_r, ~] = Optimal_2D_TrajCW(params);
        %%% Add spacecraft rotation to the optimal trajectory design
        [x_r, u_r] = Optimal_2D_TrajCW_theta(params, x_r);
        u_r = u_r/params.thrust_factor; % because |B| in Optimal_2D_TrajCW is 1 instead of thrust_factor
        %%% Save the trajectory
        interval_id = 1+N*(transferNum-1):N*transferNum;
        X_ref(:,interval_id) = x_r;
        U_ref(:,interval_id) = u_r;
    end
else
    %%% Load a previous reference trajectory with coarser time steps
    DT = 2*dt; % coarser time step
    filename = 'data/ref_traj_' + string(60*params.transfer_time) + 'min_dt=' + string(DT) + '.mat';
    while ~isfile(filename) && DT < 100*dt
        DT = DT + dt;
        filename = 'data/ref_traj_' + string(60*params.transfer_time) + 'min_dt=' + string(DT) + '.mat';
    end
    if ~isfile(filename)
        error('The desired time step is too small, start by calculating the reference trajectory for a larger dt.')
    else
        rough_traj = load(filename);
    end
    numSteps = round(DT/dt); % number of steps of size dt in each DT
    interval_id = 1:numSteps;
    %%% Linear interpolation inbetween states
    for i = 1:length(rough_traj.X_ref(1,:))-1

        %%% Define start and end point of the transfer
        x0 = rough_traj.X_ref(:, i);
        x1 = rough_traj.X_ref(:, i+1);
        for j = 1:length(x0)
            x_r(j,:) = linspace(x0(j), x1(j), numSteps+1);
        end
        X_ref(:,interval_id) = x_r(:,1:end-1);
        U_ref(:,interval_id) = rough_traj.U_ref(:, i)* ones(1, numSteps);
        interval_id = interval_id + numSteps;
    end
    x0 = x_r(:, end-1);
    x1 = rough_traj.X_ref(:, end);
    for j = 1:length(x0)
        x_f(j,:) = linspace(x0(j), x1(j), numSteps);
    end
    X_ref(:,interval_id) = x_f;
    U_ref(:,interval_id) = rough_traj.U_ref(:, end)* ones(1, numSteps)/numSteps;
end



%%% Storing data
filename = 'data/ref_traj_' + string(60*params.transfer_time) + 'min_dt=' + string(dt) + '.mat';
save(filename, 'X_ref', 'U_ref')

end