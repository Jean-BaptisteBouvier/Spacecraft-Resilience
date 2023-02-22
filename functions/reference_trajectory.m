%%% Generating the reference trajectory

%%% Based on the method from the paper:
%%% "Autonomous Optimal Trajectory Planning for Orbital Rendezvous,
%%% Satellite Inspection, and Final Approach Based on Convex Optimization"
%%% from Nicholas Ortolano, David K. Geller, and Aaron Avery.
%%% Non rotating target. Optimal trajectories for the Hill-Clohessy-Wiltshire 
%%% model in cartesian coordinates using SOCPs. 

function [X_ref, U_ref] = reference_trajectory(SOCP_Params)

dt = SOCP_Params.dt; % desired time step
N = SOCP_Params.simTimeHours*60*60/dt;
nb_transfers = length(SOCP_Params.waypoints(:,1))-1;
[n, m] = size(SOCP_Params.matrix_B);
X_ref = zeros(n, N*nb_transfers);
U_ref = zeros(m, N*nb_transfers);
SOCP_Params.numSteps = N;


for transferNum = 1:nb_transfers
    %%% Define start and end point of the transfer
    SOCP_Params.rInit = SOCP_Params.waypoints(transferNum,:)'; % x is r-bar, y is v-bar, z is cross track        
    SOCP_Params.rFinal = SOCP_Params.waypoints(transferNum+1,:)';
    %%% Calculate trajectory without the spacecraft rotation
    [~, x_r, ~] = Optimal_2D_TrajCW(SOCP_Params);
    %%% Add spacecraft rotation to the optimal trajectory design
    [x_r, u_r] = Optimal_2D_TrajCW_theta(SOCP_Params, x_r);
    u_r = u_r/SOCP_Params.thrust_factor; % because |B| in Optimal_2D_TrajCW is 1 instead of thrust_factor
    %%% Save the trajectory
    interval_id = 1+N*(transferNum-1):N*transferNum;
    X_ref(:,interval_id) = x_r;
    U_ref(:,interval_id) = u_r;
end

%%% Storing data
save('ref_traj.mat', 'X_ref', 'U_ref')

end