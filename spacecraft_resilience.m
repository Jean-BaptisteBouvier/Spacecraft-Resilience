%%% Verifies the spacecraft resilience by computing the Minkowski difference
%%% of the input sets for a chosen actuator failure.
%%% Computes also the reachable set from a given initial state.

%%% Initial dynamics are dX/dt = A*X(t) + r*R_theta*B_bar*u_bar(t)
%%% with u_bar in [0, 1]^5, R_theta a rotation matrix depending of X(t), and
%%% r is a thrust factor = thrusters_max_power/spacecraft_mass.
%%% After the loss of control over ONE thruster failure_id, the
%%% malfunctioning dynamics become:
%%% dX/dt = A*X(t) + r*R_theta*B*u(t) + r*R_theta*C*w(t) 
%%% with u in U = [0, 1]^4 and w in W = [0, 1].
%%% Following resilience theory we take the Minkowski difference of the
%%% input sets: P = BU - (-CW). The system is resilient if P is non-empty.
%%% We then calculate the reachable set of the resilient dynamics
%%% dX/dt = A*X(t) + r*R_theta*p(t)  with p(t) in P.


% Initial control matrix
B_bar = [1, 1, -1, -sqrt(2), -1;
         1, -1, -1, 0, 1];
% index of the malfunctioning thruster
failure_id = 1; % in [1,5]
% Matrix of uncontrolled thrusters
C = B_bar(:,failure_id);
% Matrix of controlled thrusters
B = B_bar; B(:,failure_id) = [];
% Matrix of the linear dynamics
A = zeros(4,4);
Omega = 0.001060922896439;
A(3,1) = 3*Omega^2; A(1,3) = 1; A(2,4) = 1;
A(3,4) = 2*Omega; A(4,3) = -2*Omega;
% Thrust factor
r = 1.5e-4; % m/s^2
% Initial state for the reachability analysis [ x, y, v_x, v_y] (m) and (m/s)
X0 = [1; 1; 0; 0]; % X0 is in the reachable set
% X0 = [100; 100; 0; 0]; % X0 is not in the reachable set
% Rotation matrix at X0
theta = atan2(X0(2), X0(1));
R_theta = [eye(2), zeros(2,2); zeros(2,2), [cos(theta), -sin(theta); sin(theta), cos(theta)]];
% Final time for reachable set computation
T = 50; % (s)
% Time step for the intermediary reachable sets
dt = 10; % (s)



%%% Calculation of inputs sets as zonotopes
BU = zonotope(sum(B'/2)', B/2); % zonotope requires the CORA and mpt toolboxes
if length(failure) > 1.5
    error('Only works for the loss of a single actuator.')
end
CW = zonotope(C/2, C/2); % If C is not a column vector, need same syntax as for BU
minus_CW = zonotope(-C/2, -C/2); % Opposite of CW
if ~BU.contains(minus_CW)
    warning('The spacecraft is not resilient to the loss of thruster no.%i.', failure_id)
end

%%% Calculation of the Minkowski difference
mD = minkDiff(BU, minus_CW); % approximated Minkowksi difference
P = reduce(mD, 'pca'); % removes empty generators
assert(isequal(P, mD)) % verification that no info is lost through the reduction
P = zonotope([zeros(2, length(P.Z(1,:))); P.Z]); % adds the first two null dimensions of the dynamics

% %%% Plotting the input sets
% figure; hold on; grid on; axis equal;
% plot(BU, [1,2], 'LineWidth', 2, 'Color', 'blue')
% plot(CW, [1,2], 'LineWidth', 2, 'Color', 'red')
% plot(P, [3, 4], 'LineWidth', 2, 'Color', 'green')
% title('Input sets BU (blue), CW (red) and their Minkowski difference P = BU - (-CW) (green)')
% set(gca,'fontsize', 18);


%%% Reachable set propagation
Phi = expm(A*dt); % transition matrix
EAB = zonotope(integral(@(t) expm(A*(dt - t))*r*R_theta*P.Z, 0, dt, 'ArrayValued', true));
Reach_Set = X0;

figure; hold on;  grid on; axis equal; 
for t = 0:dt:T
    Reach_Set = Phi*Reach_Set + EAB;
    plot(100*Reach_Set, [3, 4], 'Color', [0 t/T 1-t/T], 'LineWidth',2);
end
scatter(100*X0(3), 100*X0(4), 20, 'blue', 'filled') % factor 100 for cm/s
xlabel('v_x (cm/s)')
ylabel('v_y (cm/s)')
title('Reachable sets in the velocity space')
set(gca,'fontsize', 18);
if X0(1) == 100
    xlim([-0.3 3.5])
    ylim([-1.2 0.2])
end
if X0(1) == 1
    xlim([-2.5 1.5])
    ylim([-1 0.2])
end



% %%% Reachable set projected in the position space
% figure; hold on; grid on; axis equal;
% scatter(X0(1), X0(2), 20, 'blue', 'filled')
% plot(Reach_Set, [1, 2], 'Color', [0 t/T 1-t/T], 'LineWidth',2);
% xlabel('x (m)')
% ylabel('y (m)')
% title('Reachable set in the position space')
% set(gca,'fontsize', 18);

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other approach for resilience using only the mpt toolbox.
%%% No reachable set computation here, only the Minkowski difference.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p = length(failure_id); % number of malfunctioning thrusters
% [n, m] = size(B);
% % Undesirable thrust input set
% W = [zeros(p,1), ones(p,1)];
% % Control set
% U = [zeros(m,1), ones(m,1)];
% 
% % All possible inputs combinations
% combination_U = unique(nchoosek(U(:), m), 'rows');
% combination_W = unique(nchoosek(W(:), p), 'rows');
% 
% % Calculation all possible Bu and Cw for u in U and w in W
% CW = zeros(2^p, n);
% BU = zeros(2^(m-p), n);
% 
% id = 1;
% for i = 1:length(combination_U(:,1))
%     comb = unique(perms(combination_U(i,:)), 'rows');
%     for j = 1:length(comb(:,1))
%         u = comb(j,:);
%         Bu = B*(u');
%         BU(id,:) = Bu';
%         id = id + 1;
%     end
% end
% BU = unique(BU, 'rows');
% 
% id = 1;
% for i = 1:length(combination_W(:,1))
%     comb = unique(perms(combination_W(i,:)), 'rows');
%     for j = 1:length(comb(:,1))
%         w = comb(j,:);
%         Cw = C*(w');
%         CW(id,:) = Cw';
%         id = id + 1;
%     end
% end
% 
% %%% Transforming them into polygones
% poly_BU = Polyhedron(BU); % Polyhedron requires the mpt toolbox
% poly_CW = Polyhedron(CW);
% poly_minus_CW = -poly_CW;
% poly_P = poly_BU - poly_minus_CW;
% 
% 
% figure(1)
% hold on
% grid on
% plot(poly_BU)
% plot(poly_minus_CW)
% axis equal
% title(sprintf('Input sets BU and -CW for failure no.%i', failure_id))
% 
% 
% figure(2)
% hold on
% grid on
% plot(poly_P)
% axis equal
% title('Minkowski difference of input sets')
