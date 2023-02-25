% This function uses CVX to generate the optimal trajectory of a convex
% RPO problem. This function is limited to a non rotating target spacecraft
% and uses the CW equations in cartesian coordinates. 
%
% Author: Himmat Panag

function [eta, x, u] = Optimal_2D_TrajCW(params)
   
    N = params.numSteps; % Number of discretization steps in control input
    
    f = ones(N,1)./params.Omega^2;
    simTime = params.simTimeHours*60*60; 
    dt = simTime/N;
    A = params.matrix_A;
    B = params.matrix_B/params.thrust_factor;
    [n, m] = size(B);

    Phi = expm(A*dt);
    
    %%% Need to calculate Bd = Phi * Integral(expm(-Au)*B du from u = 0 to dt)
    % Equation 8 in 2021 Ortolano paper 
    % Method 1, Source: https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential#:~:text=is%20given%20by%20x(t,eAtx0.&text=for%20some%20T%3E0.
    % Construct matrix A2 = [-A,B;0;0]; 
    A2 = [-A,B;zeros(m,n+m)];
    A2tExp = expm(A2*dt);
    IntegralExpAdtB = A2tExp(1:n, n+1:n+m);
    Bd = Phi*IntegralExpAdtB;
    
    %%% Optimization
    cvx_begin quiet
        variable eta(N)
        variable x(n,N)
        variable u(m,N)

        minimize(f'*eta)
        subject to 
        x(:,1) == [params.rInit; params.vInit];
        x(:,N) == [params.rFinal; params.vFinal];
        
        for ii = 1:N
            norm(u(:,ii)) <= eta(ii); % thrust minimizer
            0 <= eta(ii);
            eta(ii) <= params.aMax; % max thrust
            for jj = 1:m
                0 <= u(jj,ii); % only positive thrust inputs
            end
            
            if params.coneConstraintActive 
                coneAxisUnitVec = params.coneAxisVector;
                norm(x(1:2,ii))*cos(params.coneAngle) <= coneAxisUnitVec'*x(1:2,ii);
            end 
            
            if params.usePlanarConstraints
                for kk = 1:numel(params.planarConstraint)
                    % Assert that dot(x_kk,i) > Rmin;
                    dot(x(1:2,params.planarConstraint(kk).elementNum), params.planarConstraint(kk).dirVector) >= params.planarConstraint(kk).minRadius;
                end 
            end 
        end 
        for ii = 2:N
            x(:,ii) == Phi*x(:,ii-1) + Bd*u(:,ii-1);
        end         
    cvx_end 
end 