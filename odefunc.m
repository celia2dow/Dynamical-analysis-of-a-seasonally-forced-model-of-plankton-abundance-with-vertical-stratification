function dxdt = odefunc(t, x, params)
% Edited ODE system for NPZD model
    % State variables
    N = x(1); P = x(2); Z = x(3);

    % Rename params
    a_hat = params.a_hat; 
    b = params.b; 
    c = params.c; 
    d = params.d; 
    e = params.e; 
    k_hat = params.k_hat; 
    r = params.r; 
    s_hat = params.s_hat; 
    N_0 = params.N_0; 
    alpha = params.alpha; 
    beta = params.beta; 
    gamma = params.gamma; 
    lambda = params.lambda; 
    mu = params.mu; 

    % Mixed layer depth
    if params.changing_MLD
        [M, h] = MLD_func(t);
    else
        M = params.assumed_MLD;
        h = 0;
    end

    % Vertical stratification scaling
    if params.vert_strat
        f_phiT = params.strat_func(t);
    else
        f_phiT = params.assumed_fphiT;
    end
    
    % ODE
    dNdt = -N*P*a_hat/(M*(e+N)*(b+c*P)) + ...
        r*P + beta*lambda*P^2*Z/(mu^2+P^2) + gamma*d*Z^2 + ...
        (k_hat*f_phiT + max(h,0)) * (N_0-N)/M;
    dPdt = N*P*a_hat/(M*(e+N)*(b+c*P)) - ...
        r*P - lambda*P^2*Z/(mu^2+P^2) - ...
        (s_hat*f_phiT + k_hat*f_phiT + max(h,0)) * P/M;
    dZdt = alpha*lambda*P^2*Z/(mu^2+P^2) - ...
        d*Z^2 - h*Z/M;

    dxdt = [dNdt; dPdt; dZdt];
end