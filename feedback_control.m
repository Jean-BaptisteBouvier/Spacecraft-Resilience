%%%%%% Calculation of the feedback control law of Theorem 4.
%
% This theorem only works for the loss of control authority over thruster no.4
% In this code we tune matrices K, P and Q in order to make epsilon smaller
% than rho_max - p_ref, so that P_epsilon + P_ref is a subset of P_b, the
% ball of radius rho_max.
% We also want to maximize the tracking tolerance because predictions are
% always much smaller than actual norm differences in the simulations.

clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Omega = 0.001060922896439; % (s^-1)
A = [0, 0, 1, 0; 0, 0, 0, 1; 3*Omega^2, 0, 0, 2*Omega; 0, 0, -2*Omega, 0];
B = [zeros(2,4); 1, 1, -1, -1; 1, -1, -1, 1];
C = [0; 0; -sqrt(2); 0];%/100;
r = 1.5e-4; % (m/s^2) thrust to weight ratio
L = 0.1; % Lipschitz constant of w
tau = 0.2; % (s) actuation delay
mu_A = max(eig( (A+A')/2 )); % log-norm of A
rho_max = sqrt(2) - norm(C)/sqrt(2); % radius of the ball P_b

%%%%%%%% Calculating the max thrust for the reference trajectory %%%%%%%%%%
% load('ref_traj_90min.mat')
% p_ref = 0;
% for i = 1:length(U_ref(1,:)) 
%     bu = norm(B*U_ref(:,i));
%     if bu > p_ref
%         p_ref = bu;
%     end
% end
p_ref = 4.849e-04;


%%%%%%%%%%%%%%%% Optimization over feedback gain K %%%%%%%%%%%%%%%%%%%%%
% for tau = 0.2s
K = 472*[1,1,1,1; 1,-1,1,-1; -1,-1,-1,-1; -1,1,-1,1]; % tol = 0.154 mm  
A_tilde = A - r*B*K;

%%%%%%%%%%% Optimization over Lyapunov matrices (P, Q) %%%%%%%%%%%%%%%%
% Scaling Q has no influence.
% Adding non-diagonal terms to Q, or modifying the diagonal terms of Q
% reduces the tracking tolerance for a same epsilon.
% Best choice is Q = I. 
% P cannot be modified once Q is chosen because they verify the Lyapunov
% equation.
Q = eye(4); 
P = lyap(A_tilde', Q);


%%%%%%%%%%%%% Calculation of the feedback parameters %%%%%%%%%%%%%%%%%%
alpha = min(eig(Q))/(2*max(eig(P)));
beta = r*sqrt(max(eig(P)))*norm(C)*L*tau;
gamma = r*norm(B*K)*(exp(mu_A*tau)-1)/mu_A;
eps = norm(B*K)*beta*(1+gamma)/(alpha*sqrt(min(eig(P)))) + gamma*norm(C)*L*tau;
tol = beta*(1+gamma)/(alpha*sqrt(min(eig(P)))); % tracking tolerance


disp('epsilon = ' + string(eps) + '    rho_max - p_ref = ' + string(rho_max - p_ref) + '    tracking tol = ' + string(tol))


