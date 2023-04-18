%%% Function to generating the undesirable thrust input w
%%% Inputs:
%%% is_Lipschitz: (bool) [true] w is Lipschitz  [false] w is bang-bang
%%% w_max: (real in [0,1]) saturation limit for the undesirable input
%%% L: (real) Lipschitz constant, only used if is_Lipschitz
%%% num_bang_per_hour: (integer) average numbers of bangs per hour, only used if ~is_Lipschitz
%%% params: structure containing all the parameters of the simulation

%%% If W is not saved in the data folder it is created with unit amplitude
%%% and stored in the data folder.
%%% Then W is scaled by w_max for use.

function W = undesirable_input(is_Lipschitz_w, w_max, L, num_bang_per_hour, params)

dt = params.dt;
N = params.transfer_time*60*60/dt;
p = length(params.matrix_C(1,:)); % number of malfunctioning thrusters
nb_transfers = length(params.waypoints(:,1))-1;

data_path = 'data/'; % path to the data folder from the main


if is_Lipschitz_w 
    
	W_Lip_filename = 'W_'+ string(60*params.transfer_time) + 'min_Lip=' + string(L) + '_dt=' +string(dt)+ '.mat';
	disp('w has a Lipschitz constant L = ' + string(L) + ' and w_max = ' + string(w_max) + 'm/s^2')
    
	if isfile(data_path+W_Lip_filename) % signal already exists
		load(data_path+W_Lip_filename, 'W');
        
    else % build signal with unit amplitude
        W = zeros(p, N*nb_transfers);
		W(:,1) = zeros(p,1);
		for i = 2:N*nb_transfers
			W(:,i) = W(:,i-1) + dt*L*2*(rand(p,1)-0.5);
			W(:,i) = ((W(:,i) <= 1).*(W(:,i) >= 0)).*W(:,i) + (W(:,i) > 1); % clamp in [0,1]
		end
		save(data_path+W_Lip_filename, 'W');
	end
else % w is bang-bang
    if round(num_bang_per_hour) ~= num_bang_per_hour
        warning('The number of bangs per hour will be rounded')
        num_bang_per_hour = round(num_bang_per_hour);
    end
    
	W_bang_filename = 'W_' + string(60*params.transfer_time) + 'min_bang=' + string(num_bang_per_hour)+ 'per_hour_dt=' + string(dt) + '.mat';
	disp('w is bang-bang with an average of '+string(num_bang_per_hour)+' bangs per hour and w_max = ' + string(w_max) + 'm/s^2')
    
	if isfile(data_path+W_bang_filename) % signal already exists
		load(data_path+W_bang_filename, 'W');
        
    else % build signal with unit amplitude
        W = zeros(p, N*nb_transfers);
        bang_probability = num_bang_per_hour * dt/3600;
		W(:,1) = ones(p,1);
		for i = 2:N*nb_transfers
			if rand > 1-bang_probability
				W(:,i) = W(:,i-1) == zeros(p,1);
			else
				W(:,i) = W(:,i-1);
			end
		end
		save(data_path+W_bang_filename, 'W');
	end
end

W = W*w_max;

end