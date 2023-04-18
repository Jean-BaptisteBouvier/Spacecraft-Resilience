%%% Function to track a given trajectory X_ref of dynamics 
%%% d/dt X_ref(t) = A X_ref(t) + B U_ref(t),  t >= 0. 
%%% The tracker is affected by input delay and perturbations.
%%% The other inputs are:
%%% time series of the undesirable inputs W,
%%% constant known time delay h and,
%%% state starts at X0, which is equal to X_ref only for the first transfer

function [X, U] = Lechappe_tracking(X0, X_ref, U_ref, W, h, params)

n = length(X_ref(:,1)); % nb of states
m = length(U_ref(:,1)); % nb of control inputs
dt = params.dt; % time step
N = params.transfer_time*60*60/dt; % nb of steps
delta_id = round(h/dt); % nb of time steps corresponding to the actuation delay

%%% Creating the ODE
A = params.matrix_A;
B = params.matrix_B;
C = params.matrix_C;

%%% Storing the tracking data
X = zeros(n, N); X(:,1) = X0;
U = zeros(m, N);

%%% Constant gain matrix obtained from `feedback_control.m`
K = 472*[1,1,1,1; 1,-1,1,-1; -1,-1,-1,-1; -1,1,-1,1];

%%% Parameters for the linear optimization
scaling = 1e6; % to help linprog converge
options = optimoptions('linprog', 'Display', 'off');
f = ones(m,1);

for i = 1:N-1
    t = i*dt;
    
    if t <= h + dt*1e-5 % not enough time steps for the Lechappe predictor
        U(:,i) = U_ref(:,i);
    else
        %%% Calculation of the integral term of the Lechappe predictor
        intg = 0;
        for j = 0:delta_id-1
            tau = (t - 2*h) + j*dt;
            id = i - delta_id + j;
            theta = atan2(X(2,id), X(1,id));
            R_theta = [zeros(2,4); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
            intg = intg + dt*expm(A*(t - h - tau))*R_theta*(B*U(:,id) + C*W(:,id));
        end
        X_prediction = expm(A*h)*X(:,i-delta_id) + intg; % Lechappe predictor
        
        theta = atan2(X(2,i), X(1,i));
        R_theta_inv = [eye(2), zeros(2,2); zeros(2,2), [cos(theta), sin(theta); -sin(theta), cos(theta)]];
        
        %%% Calculation of control input u_w counteracting Cw
        u_w = linprog(f, [], [], B*scaling, -C*W(:,i-delta_id)*scaling, zeros(m,1), ones(m,1), options);
        if isempty(u_w)
            u_w = zeros(4,1);
        end
        %%% Calculation of the feedback control input u_eps
        u_eps = linprog(f, [], [], B*scaling, R_theta_inv*B*K*(X_ref(:,i) - X_prediction)*scaling, zeros(m,1), [], options);
        if isempty(u_eps)
            u_eps = zeros(4,1);
        end
        %%% Total control input
        U(:,i) = U_ref(:,i) + u_eps + u_w; 
         
        %%% Input bounds [0, 1]
        U(:,i) = (U(:,i) >= 0).*(U(:,i) <= 1).*U(:,i) + (U(:,i) > 1);        
    end
   
    %%% Propagation of the state
    theta = atan2(X(2,i), X(1,i));
    R_theta = [zeros(2,4); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
    X(:,i+1) = X(:,i) + dt*(A*X(:,i) + R_theta*B*U(:,i) + R_theta*C*W(:,i));
end


end