%%% Function to track a given trajectory X_ref of dynamics 
%%% d/dt X_ref = A X_ref + B U_ref 
%%% The tracker is affected by input delay and perturbations

%%% Time series of the undesirable inputs W
%%% Constant known time delay is h
%%% State starts at X0, which is equal to X_ref only for the first transfer

function [X_Lechappe, U_Lechappe] = Lechappe_tracking(X0, X_ref, U_ref, W, h, params)

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
end


end




%%%%%%%%%% Previous version

% function [X_Lechappe, U_Lechappe] = Lechappe_tracking(X0, X_ref, U_ref, W, h, params)
% 
% n = length(X_ref(:,1)); % nb of states
% m = length(U_ref(:,1)); % nb of control inputs
% dt = params.dt;
% N = params.simTimeHours*60*60/dt;
% 
% %%% Creating the ODE
% A = params.matrix_A;
% B = params.matrix_B;
% C = params.matrix_C;
% 
% %%% Tracking
% X_Lechappe = zeros(n, N); X_Lechappe(:,1) = X0;
% U_Lechappe = zeros(m, N);
% delta_id = round(h/dt);
% X_goal = [params.rFinal; params.vFinal]; 
% eps = 10e-5; % 10 (cm)
% stop_transfer = false;
% 
% %%% Parameters for the linear optimization
% scaling = 1e6; % to help linprog converge
% options = optimoptions('linprog', 'Display', 'off');
% f = ones(m,1);
% 
% for i = 1:N-1
%     t = i*dt;
%     
%     if t <= h % not enough time steps for the Lechappe predictor
%         U_Lechappe(:,i) = 0;
%     else
%         intg = 0; % integral term in the Lechappe predictor
%         for j = 0:delta_id-1
%             tau = (t - 2*h) + j*dt;
%             id = i - delta_id + j;
%             theta = atan2(X_Lechappe(2,id), X_Lechappe(1,id));
%             R_theta = [zeros(2,4); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
%             intg = intg + dt*expm(A*(t - h - tau))*R_theta*(B*U_Lechappe(:,id) + C*W(:,id));
%         end
%         X_prediction = expm(A*h)*X_Lechappe(:,i-delta_id) + intg;
%         u_w = linprog(f, [], [], B*scaling, -C*W(:,i-delta_id)*scaling, zeros(m,1), ones(m,1), options);
%         theta = atan2(X_Lechappe(2,i), X_Lechappe(1,i));
%         
%         %%% Rotating gain
% %         R_theta = [eye(2), zeros(2,2); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
% %         K = lqr(A, R_theta*B, eye(n), eye(m));
% %         if stop_transfer
% %             U_Lechappe(:,i) = U_ref(:,i - delta_id) + K*(X_goal - X_prediction) + u_w;
% %         else
% %             U_Lechappe(:,i) = U_ref(:,i - delta_id) + K*(X_ref(:,i) - X_prediction) + u_w;
% %         end
%         
%         %%% Constant gain obtained from `feedback_control.m`
%         R_theta_inv = [eye(2), zeros(2,2); zeros(2,2), [cos(theta), sin(theta); -sin(theta), cos(theta)]];
%         K = 500*[1,1,1,1; 1,-1,1,-1; -1,-1,-1,-1; -1,1,-1,1];
%         
% %         if stop_transfer
% %             u_eps = linprog(f, [], [], B*scaling, R_theta_inv*B*K*(X_goal - X_prediction)*scaling, zeros(m,1), ones(m,1), options);
% %         else
%             u_eps = linprog(f, [], [], B*scaling, R_theta_inv*B*K*(X_ref(:,i) - X_prediction)*scaling, zeros(m,1), ones(m,1), options);     
% %         end
%         if isempty(u_eps) % no solution for u_eps to follow X_ref, instead try to go to X_goal
%             u_eps = linprog(f, [], [], B*scaling, R_theta_inv*B*K*(X_goal - X_prediction)*scaling, zeros(m,1), ones(m,1), options);
%         elseif isempty(u_w)
%             u_w = zeros(4,1);
%         end
%         U_Lechappe(:,i) = U_ref(:,i - delta_id) + u_eps + u_w; 
%         
%         
%         %%% Input bounds [0, 1]
%         U_Lechappe(:,i) = (U_Lechappe(:,i) >= 0).*(U_Lechappe(:,i) <= 1).*U_Lechappe(:,i) + (U_Lechappe(:,i) > 1);        
%     end
%    
%     theta = atan2(X_Lechappe(2,i), X_Lechappe(1,i));
%     R_theta = [zeros(2,4); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
%     X_Lechappe(:,i+1) = X_Lechappe(:,i) + dt*(A*X_Lechappe(:,i) + R_theta*B*U_Lechappe(:,i) + R_theta*C*W(:,i));
% %     scatter(X_Lechappe(1,i+1)*1e3, X_Lechappe(2,i+1)*1e3, 20, 'red', 'filled')
% 
%     if ~stop_transfer && (norm( X_Lechappe(:,i+1) - X_goal ) < eps || norm( X_ref(:,i+1) - X_goal ) < eps) % if close enough, switch from trajectory tracking to CV to X_goal
%         stop_transfer = true;
%     end
% end
% 
% 
% end