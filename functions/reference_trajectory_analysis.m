%%% Analysis of the reference trajectory (X_ref, U_ref)

function reference_trajectory_analysis(params, X_ref, U_ref)

nb_transfers = length(params.waypoints(:,1))-1;
N = length(X_ref(1,:));
dt = params.dt;
time = [1:N]*dt/3600; % [hours]
x_limit = nb_transfers*params.simTimeHours;




%%% Plotting the trajectories, waypoints, KOS   with X and Y inverted
R_KOS = 50; % [m] radius of the Keep-Out-Sphere from the Restore-L mission
colors = get(gca,'colororder');
hold on; grid on;
plot( R_KOS*cos(0:0.01:2*pi), R_KOS*sin(0:0.01:2*pi), 'Color', colors(3,:), 'LineWidth', 2) 
scatter(0,0,50,'red','filled') % target spacecraft
plot(X_ref(2,:)*1e3, X_ref(1,:)*1e3, 'Color', colors(1,:), 'LineWidth', 2); % plotting reference trajectory
for transferNum = 1:nb_transfers % plotting waypoints
    scatter(params.waypoints(transferNum,2)*1e3, params.waypoints(transferNum,1)*1e3, 50, colors(5,:), 'filled')
end
set(gca,'fontsize', 18);
axis equal
xlabel('Y (m)'); ylabel('X (m)');
xlim([-420 420]); ylim([-150 150]); % Cropping the trajectory





%%% Reference thrust
B = params.matrix_B / params.thrust_factor;
p_ref_max = 0;
for i = 1:N 
    p_ref = norm(B*U_ref(:,i));
    if p_ref > p_ref_max
        p_ref_max = p_ref;
    end
end
disp('The maximal input norm for the reference trajectory is  max ||p_ref|| = ' + string(p_ref_max))

figure
hold on
grid on
for i = 1:10:N % too many datapoints for scatter
    p_ref = B*U_ref(:,i);
    scatter(p_ref(3), p_ref(4), 30, 'blue', 'filled')
end
xlabel('p_{ref}(3)   (m/s^2)')
ylabel('p_{ref}(4)   (m/s^2)')
set(gca,'fontsize', 18);



%%% Orientation on the reference orbit
figure
hold on
grid on
for i = 1:100:N % too many datapoints for scatter
    theta = atan2(X_ref(2,i), X_ref(1,i))*180/pi; % rad to deg
    if theta < 0
        theta = theta + 360;
    end
    scatter(time(i), theta, 20, colors(1,:), 'filled')
end
for transferNum = 1:nb_transfers+1 % plotting waypoints
    theta = atan2(params.waypoints(transferNum,2), params.waypoints(transferNum,1))*180/pi; % rad to deg
    if theta < 0
        theta = theta + 360;
    end
    scatter(params.simTimeHours*(transferNum-1), theta, 40, colors(5,:), 'filled')
    if theta == 0
        scatter(params.simTimeHours*(transferNum-1), 360, 40, colors(5,:), 'filled')
    end
end
xlim([0 x_limit])
ylim([0 360])
set(gca,'fontsize', 18);
xlabel('time (hours)')
ylabel('\theta (deg)')
yticks([0, 90, 180, 270, 360])
yticklabels({'0^\circ','90^\circ','180^\circ','270^\circ','360^\circ'})




%%% Reference velocity
v_ref = zeros(1, N);
for i = 1:N
    v_ref(i) = norm(X_ref(3:4, i));
end

figure
hold on
grid on
plot(time, v_ref*1e5, 'LineWidth', 2)
xlabel('time (hours)')
ylabel('velocity (cm/s)')
xlim([0 x_limit])
set(gca,'fontsize', 18);





%%% Reference fuel consumption
V_exit = params.V_exit; % [m/s] exit velocity of ions in PPS-1350
M_ref = zeros(1, N);
M_ref(1) = params.mass;

for i = 2:N
    M_ref(i) = M_ref(i-1) - dt*M_ref(i-1)*norm(U_ref(:,i))/V_exit;
end
% Get mass of fuel consumed instead of spacecraft mass
M_ref = params.mass - M_ref;
disp('m_ref = ' + string(M_ref(end)) + 'kg.')

figure
hold on
grid on
ylabel('Mass of fuel consumed (kg)')
xlabel('time (hours)')
plot(time, M_ref, 'LineWidth', 2)
xlim([0 x_limit])
set(gca,'fontsize', 18);


%%% Norm of the control inputs
ref_control = zeros(1, N);
for i = 1:N
    ref_control(i) = norm(U_ref(:,i));
end

figure
hold on
grid on
plot(time, U_ref(1,:), 'LineWidth', 2)
plot(time, U_ref(2,:), 'LineWidth', 2)
plot(time, U_ref(3,:), 'LineWidth', 2)
plot(time, U_ref(4,:), 'LineWidth', 2)
thrusters = 1:5; thrusters(thrusters == params.failure) = [];
legend('u_{ref}^'+string(thrusters(1)), 'u_{ref}^'+string(thrusters(2)), 'u_{ref}^'+string(thrusters(3)), 'u_{ref}^'+string(thrusters(4)))
xlabel('time (hours)')
ylabel('input (m/s^2)')
xlim([0 x_limit])
set(gca,'fontsize', 18);


%%% Magnitude of the reference control
figure
hold on
grid on
plot(time, ref_control, 'LineWidth', 2)
xlabel('time (hours)')
ylabel('input (m/s^2)')
xlim([0 x_limit])
set(gca,'fontsize', 18);



end