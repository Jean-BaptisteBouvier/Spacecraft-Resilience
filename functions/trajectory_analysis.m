%%% Analysis of the trajectory (X, U) with respect to the reference
%%% trajectory (X_ref, U_ref) with undesirable thrust W

function trajectory_analysis(params, X_ref, U_ref, W, X, U)

nb_transfers = length(params.waypoints(:,1))-1;
N = length(X_ref(1,:));
dt = params.dt;
time = [1:N]*dt/3600; % [hours]
x_limit = nb_transfers*params.simTimeHours;




%%% Plotting the trajectories, waypoints, KOS   with X and Y inverted
R_KOS = 50; % [m] radius of the Keep-Out-Sphere from the Restore-L mission
figure
colors = get(gca,'colororder');
hold on; grid on;
plot( R_KOS*cos(0:0.01:2*pi), R_KOS*sin(0:0.01:2*pi), 'Color', colors(3,:), 'LineWidth', 2) 
scatter(0,0,50,'red','filled') % target spacecraft
plot(X_ref(2,:)*1e3, X_ref(1,:)*1e3, 'Color', colors(1,:), 'LineWidth', 2); % plotting reference trajectory
plot(X(2,:)*1e3, X(1,:)*1e3, 'Color', colors(2,:), 'LineWidth', 2); % plotting tracking trajectory
for transferNum = 1:nb_transfers % plotting waypoints
    scatter(params.waypoints(transferNum,2)*1e3, params.waypoints(transferNum,1)*1e3, 50, colors(5,:), 'filled')
end
set(gca,'fontsize', 18);
axis equal
xlabel('Y (m)'); ylabel('X (m)');
xlim([-420 420]); ylim([-150 150]); % Cropping the trajectory





%%% Speed
[speed_ref, speed_Lechappe] = deal(zeros(1, N));
for i = 1:N
    speed_ref(i) = norm(X_ref(3:4, i));
    speed_Lechappe(i) = norm(X(3:4, i));
end

figure
hold on
grid on
plot(time, speed_ref*1e5, 'LineWidth', 2)
plot(time, speed_Lechappe*1e5, 'LineWidth', 2)
xlabel('time (hours)')
ylabel('velocity (cm/s)')
legend('v_{ref}', 'v')
xlim([0 x_limit])
set(gca,'fontsize', 18);


%%% Norm error
norm_dif = zeros(1, N);
step_around = 5; % number of steps to look around for minimal position error

for i = 1:N
    x_ref = X_ref(:,i);
    min_dif = norm(x_ref - X(:,i) );
    if i > step_around
        for j = i-1:-1:i-step_around
            dif = norm(x_ref - X(:,j) );
            if dif < min_dif
                min_dif = dif;
            end
        end
    end
    if N - i > step_around
        for j = i+1:i+step_around
            dif = norm(x_ref - X(:,j) );
            if dif < min_dif
                min_dif = dif;
            end
        end
    end
    norm_dif(i) = min_dif*1e6; % [mm]

end
disp('Average norm error ' + string(mean(norm_dif)) + '  maximal error ' + string(max(norm_dif)))

    
    
%%% Position Error
pos_dif = zeros(1, N);
step_around = 5; % number of steps to look around for minimal position error

for i = 1:N
    x_ref = X_ref(1:2, i);
    min_dif = norm(x_ref - X(1:2, i) );
    if i > step_around
        for j = i-1:-1:i-step_around
            dif = norm(x_ref - X(1:2, j) );
            if dif < min_dif
                min_dif = dif;
            end
        end
    end
    if N - i > step_around
        for j = i+1:i+step_around
            dif = norm(x_ref - X(1:2, j) );
            if dif < min_dif
                min_dif = dif;
            end
        end
    end
    pos_dif(i) = min_dif*1e6; % [mm]
end
disp('Average position error ' + string(mean(pos_dif)) + 'mm  maximal error ' + string(max(pos_dif)) + 'mm')

figure
hold on
grid on
plot(time, pos_dif, 'LineWidth', 2)
xlabel('time (hours)')
ylabel('position error (mm)')
xlim([0 x_limit])
set(gca,'fontsize', 18);


%%% Fuel consumption
V_exit = params.V_exit; % [m/s] exit velocity of ions in PPS-1350
[M_ref, M_u, M_w] = deal(zeros(1, N));
[M_ref(1), M_u(1), M_w(1)] = deal(params.mass);

for i = 2:N
    M_ref(i) = M_ref(i-1) - dt*M_ref(i-1)*norm(U_ref(:,i))/V_exit;
    M_u(i) = M_u(i-1) - dt*M_u(i-1)*norm(U(:,i))/V_exit;
    M_w(i) = M_w(i-1) - dt*M_w(i-1)*norm(W(:,i))/V_exit;
end
% Get mass of fuel consumed instead of spacecraft mass
M_ref = params.mass - M_ref;
M_u = params.mass - M_u;
M_w = params.mass - M_w;
disp('m_ref = ' + string(M_ref(end)) + 'kg    m_u = ' + string(M_u(end)) + 'kg   m_w = ' + string(M_w(end)) + 'kg.')

figure
hold on
grid on
ylabel('Mass of fuel consumed (kg)')
xlabel('time (hours)')
plot(time, M_ref, 'LineWidth', 2)
plot(time, M_u, 'LineWidth', 2)
plot(time, M_w, 'LineWidth', 2)
legend('m_{ref}', 'm_u', 'm_w')
xlim([0 x_limit])
set(gca,'fontsize', 18);


%%% Control
[ref_control, Lechappe_control, W_control] = deal(zeros(1, N));
for i = 1:N
    ref_control(i) = norm(U_ref(:,i));
    Lechappe_control(i) = norm(U(:,i));
    W_control(i) = norm(W(:,i));
end

figure
hold on
grid on
plot(time, ref_control, 'LineWidth', 2)
plot(time, Lechappe_control, 'LineWidth', 2)
plot(time, W_control, 'LineWidth', 2)
legend('||u_{ref}||', '||u||', '||w||')
xlabel('time (hours)')
ylabel('input (m/s^2)')
xlim([0 x_limit])
set(gca,'fontsize', 18);


figure
hold on
grid on
plot(time, U(1,:), 'LineWidth', 2)
plot(time, U(2,:), 'LineWidth', 2)
plot(time, U(3,:), 'LineWidth', 2)
plot(time, U(4,:), 'LineWidth', 2)
thrusters = 1:5; thrusters(thrusters == params.failure) = [];
legend('u_'+string(thrusters(1)), 'u_'+string(thrusters(2)), 'u_'+string(thrusters(3)), 'u_'+string(thrusters(4)))
xlabel('time (hours)')
ylabel('input (m/s^2)')
xlim([0 x_limit])
set(gca,'fontsize', 18);




end