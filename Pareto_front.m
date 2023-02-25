%%% Pareto front
%%%
%%% Authors: Jean-Baptiste Bouvier and Himmat Panag.


clearvars
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%% User inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_dist = 0.8; % [m] maximal distance allowed between trajectory and reference

%%% Optimization variables
delay = 0:10; % [s] range of actuation delays to test
w_accuracy = 0.005; % accuracy on the determination of the maximal saturation value for w where the tracking converges


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 0.1; % Lipschitz constant of w
num_bang_per_hour = 10; % average number of bangs per hour when w is bang-bang
dt = 1; % [s] time step must be smaller than the delay
mass = 600; % [kg] spacecraft mass
thrust = 90e-3; % [N] max thrust of each thruster (90 mN for the PPS-1350)
failure = 4; % id of the malfunctioning thruster in {1,2,3,4,5}.  Only resilient to no. 4
transfer_time = 1.5; % [hours] transfer time between each waypoints
waypoints = [0,200; 0,80; 80,0; 0,-80; -80,0; 0,80].*1e-3; % [km] waypoints positions


%%%%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Pareto front calculation')
disp('Loss of control authority over thruster no. ' + string(failure)+ '   Transfer time = ' + string(60*transfer_time) + 'min')


%%% Add path to the folders of functions and data
addpath 'functions'
addpath 'data'
%%% Setting up all the parameters
params = parameters_setup(dt, failure, mass, thrust, transfer_time, waypoints);

%%% Generating the reference trajectory if not already done
filename = 'data/ref_traj_' + string(60*transfer_time) + 'min.mat';
if isfile(filename)
    load(filename)
else
    [X_ref, U_ref] = reference_trajectory(params);
end


%%% Calculation of the Pareto front 
disp('w has a Lipschitz constant L = ' + string(L))
Lipschitz_front = dichotomy(delay, true, w_accuracy, L, num_bang_per_hour, max_dist, params, X_ref, U_ref);
disp('w is bang-bang and with an average of ' + string(num_bang_per_hour) + ' bangs per hour.')
bang_front = dichotomy(delay, false, w_accuracy, L, num_bang_per_hour, max_dist, params, X_ref, U_ref);
save('data/pareto_front.mat', 'bang_front', 'Lip_front')


%%% Plotting Pareto front
figure; hold on; grid on;
for i = 1:length(delay)
    scatter(delay(i), Lipschitz_front(i), 30, 'blue', 'filled')
    scatter(delay(i), bang_front(i), 30, 'red', 'filled')
end
legend('Lipschitz', 'bang-bang')
xlabel('actuation delay \tau (s)')
ylabel('saturation of w')
set(gca,'fontsize', 18);








%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Calculation of the Pareto front using the
%%% dichotomy approach on the saturation value of w 
function front = dichotomy(delay, is_Lipschitz_w, w_accuracy, L, num_bang_per_hour, max_dist, params, X_ref, U_ref)


%%% Store the maximal w_max for each delay at which the tracking never
%%% exceeds max_dist.
front = zeros(1, length(delay));

for id_delay = 1:length(delay)
    
    tau = delay(id_delay);
    disp('Actuation delay = '+ string(tau) + 's')
    w_low = 0; w_up = 1;
    if id_delay > 1
        w_up = front(id_delay-1); % at a higher delay cannot overcome higher w
    end
        
    while w_up - w_low > w_accuracy
        w_max = (w_up + w_low)/2;
        %%% Generating the undesirable inputs
        W = undesirable_input(is_Lipschitz_w, w_max, L, num_bang_per_hour, params);
        %%% Verifying whether the tracking diverges
        diverge = tracking_loop(X_ref, U_ref, W, tau, params, max_dist);
       
        if diverge
            disp('    For w_max = ' + string(w_max) + 'm/s^2, tracking diverged.')
            w_up = w_max;
        else
            disp('    For w_max = ' + string(w_max) + 'm/s^2, tracking converged.')
            w_low = w_max;
        end
    end
    front(id_delay) = (w_up + w_low)/2;   
end

end




%%% Verification of whether the controller can track the reference
%%% trajectory with the given parameters: W and tau.
%%% returns boolean diverge
function diverge = tracking_loop(X_ref, U_ref, W, tau, params, max_dist)

N = params.simTimeHours*60*60/params.dt;
diverge = false;
for transferNum = 1:2 % if the first two transfers are successful, so are the following ones by symmetry
    interval_id = 1+N*(transferNum-1):N*transferNum;
    %%% Selecting the transfer start point
    if transferNum == 1
        x0 = X_ref(:,1);
    else
        x0 = x;
    end    
    %%% Trajectory tracking
    [x, diverge] = tracking(x0, X_ref(:, interval_id), U_ref(:, interval_id), W(:, interval_id), tau, params, max_dist);
    if diverge
        break
    end
end

end



%%% Generating the undesirable inputs
function W = undesirable_input(is_Lipschitz_w, w_max, L, num_bang_per_hour, params)

dt = params.dt;
N = params.simTimeHours*60*60/dt;
p = length(params.matrix_C(1,:)); % number of malfunctioning thrusters
nb_transfers = length(params.waypoints(:,1))-1;

W = zeros(p, N*nb_transfers);
if is_Lipschitz_w
    W(:,1) = w_max*ones(p,1);
    for i = 2:N*nb_transfers
        W(:,i) = W(:,i-1) + dt*L*2*(rand(p,1)-0.5);
        W(:,i) = ((W(:,i) <= w_max).*(W(:,i) >= 0)).*W(:,i) + (W(:,i) > w_max)*w_max;
    end
else % w is bang-bang
    bang_probability = num_bang_per_hour * dt/3600;
    W(:,1) = w_max*ones(p,1);
    for i = 2:N*nb_transfers
        if rand > 1-bang_probability
            W(:,i) = (W(:,i-1) == zeros(p,1))*w_max;
        else
            W(:,i) = W(:,i-1);
        end
    end
end

end


%%% Modified Lechappe tracking
%%% Faster code for the Pareto front calculation without data storage and
%%% stopping criterion: ||X_ref(t) - X(t)|| >= max_dist
%%% returns last state x and boolean diverge.

function [x, diverge] = tracking(X0, X_ref, U_ref, W, h, params, max_dist)

n = length(X_ref(:,1)); % nb of states
m = length(U_ref(:,1)); % nb of control inputs
dt = params.dt;
N = params.simTimeHours*60*60/dt;

%%% Creating the ODE
A = params.matrix_A;
B = params.matrix_B;
C = params.matrix_C;

%%% Tracking
X_Lechappe = zeros(n, N); X_Lechappe(:,1) = X0;
U_Lechappe = zeros(m, N);
delta_id = round(h/dt);
diverge = false;

%%% Constant gain obtained from `feedback_control.m`
K = 472*[1,1,1,1; 1,-1,1,-1; -1,-1,-1,-1; -1,1,-1,1];

%%% Parameters for the linear optimization
scaling = 1e6; % to help linprog converge
options = optimoptions('linprog', 'Display', 'off');
f = ones(m,1);

for i = 1:N-1
    t = i*dt;
    
    if t <= h % not enough time steps for the Lechappe predictor
        U_Lechappe(:,i) = 0;
    else
        intg = 0; % integral term in the Lechappe predictor
        for j = 0:delta_id-1
            tau = (t - 2*h) + j*dt;
            id = i - delta_id + j;
            theta = atan2(X_Lechappe(2,id), X_Lechappe(1,id));
            R_theta = [zeros(2,4); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
            intg = intg + dt*expm(A*(t - h - tau))*R_theta*(B*U_Lechappe(:,id) + C*W(:,id));
        end
        X_prediction = expm(A*h)*X_Lechappe(:,i-delta_id) + intg;
        
        theta = atan2(X_Lechappe(2,i), X_Lechappe(1,i));
        R_theta_inv = [eye(2), zeros(2,2); zeros(2,2), [cos(theta), sin(theta); -sin(theta), cos(theta)]];
        
        %%% Double linear optimization because u_eps might not work
        u_w = linprog(f, [], [], B*scaling, -C*W(:,i-delta_id)*scaling, zeros(m,1), ones(m,1), options);
        u_eps = linprog(f, [], [], B*scaling, R_theta_inv*B*K*(X_ref(:,i) - X_prediction)*scaling, zeros(m,1), ones(m,1), options);
        
        if isempty(u_eps)
            u_eps = zeros(4,1);
        elseif isempty(u_w)
            u_w = zeros(4,1);
        end
        U_Lechappe(:,i) = U_ref(:,i - delta_id) + u_eps + u_w; 
        
        %%% Input bounds [0, 1]
        U_Lechappe(:,i) = (U_Lechappe(:,i) >= 0).*(U_Lechappe(:,i) <= 1).*U_Lechappe(:,i) + (U_Lechappe(:,i) > 1);        
    end
   
    theta = atan2(X_Lechappe(2,i), X_Lechappe(1,i));
    R_theta = [zeros(2,4); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
    X_Lechappe(:,i+1) = X_Lechappe(:,i) + dt*(A*X_Lechappe(:,i) + R_theta*B*U_Lechappe(:,i) + R_theta*C*W(:,i));
    
    if norm(X_Lechappe(1:2,i+1) - X_ref(1:2,i+1)) > max_dist*1e-3 % all in km
        diverge = true;
        break
    end
end
x = X_Lechappe(:,i+1);

end


