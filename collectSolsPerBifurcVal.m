function Xstore_wStability = collectSolsPerBifurcVal(Xstore, param_i, tol, vec, Xstore_wStability)
% Function for collecting all stable and unstable steady state and limit 
% cycle solutions and storing them in the same structure with inputs:
%   Xstore (the long-term solutions for compartment X given an array of bifurcation parameter values),
%   param_i (the index for the bifurcation parameter value),
%   tol ( tolerance for the solution being steady in the long-term),
%   vec (array of parameter values for bifurcation parameter),
%   and Xstore_wStability (structure storing unique, ordered values and stabilities of long-term solutions for fsolve and ode45 solutions).
% and output:
%   Xstore_wStability.

    % FSOLVE solution
    ssXics1 = ~isnan(Xstore(param_i,:,1));                  % Initial conditions for which FSOLVE solution exists
    Xstore_1 = squeeze(Xstore(param_i,ssXics1,1:2));        % Extract columns for given bifurcation parameter pertaining to the FSOLVE solution
    if sum(ssXics1) == 1
        Xstore_1 = Xstore_1';
    end
    [Xstore_stable_1, Xstore_unstable_1] = findStabUnstabSols(Xstore_1, tol);

    % ODE45 solution: steady states
    ssXics2 = isnan(Xstore(param_i,:,4)) & ...              % Initial conditions for which ODE45 solution is a steady state
        ~isnan(Xstore(param_i,:,5));                        % and for which the limit exists
    Xstore_ss2 = squeeze(Xstore(param_i,ssXics2,[3,6]));    % Extract columns for given bifurcation parameter pertaining to ODE45 derived steady state solutions
    if sum(ssXics2) == 1
        Xstore_ss2 = Xstore_ss2';
    end
    [Xstore_stable_ss2, Xstore_unstable_ss2] = findStabUnstabSols(Xstore_ss2, tol);
    
    % ODE45 solution: limit cycles
    lcXics2 = ~isnan(Xstore(param_i,:,4)) & ...             % Initial conditions for which ODE45 solution is oscillatory
        ~isnan(Xstore(param_i,:,5));                        % and for which the limit exists
    Xstore_lc2 = squeeze(Xstore(param_i,lcXics2,[3,4,6]));  % Extract columns for given bifurcation parameter pertaining to ODE45 derived limit cycle solutions
    if sum(lcXics2) == 1
        Xstore_lc2 = Xstore_lc2';
    end
    [Xstore_stable_lc2, Xstore_unstable_lc2] = findStabUnstabSols(Xstore_lc2, tol);
    
    % COLLATE ALL SOLUTIONS
    Xstore_wStability.stable.sol_1 = ...   % fsolve steady state solutions 
        [Xstore_wStability.stable.sol_1; ...
            [vec(param_i) * ones(size(Xstore_stable_1,1),1), Xstore_stable_1] ...
        ];  
    Xstore_wStability.unstable.sol_1 = ...  
        [Xstore_wStability.unstable.sol_1; ...
            [vec(param_i) * ones(size(Xstore_unstable_1,1),1), Xstore_unstable_1] ...
        ]; 
    Xstore_wStability.stable.sol_ss2 = ... % ode45 steady state solutions 
        [Xstore_wStability.stable.sol_ss2; ...
            [vec(param_i) * ones(size(Xstore_stable_ss2,1),1), Xstore_stable_ss2] ...
        ];   
    Xstore_wStability.unstable.sol_ss2 = ... 
        [Xstore_wStability.unstable.sol_ss2; ...
            [vec(param_i) * ones(size(Xstore_unstable_ss2,1),1), Xstore_unstable_ss2] ...
        ];
    Xstore_wStability.stable.sol_lc2 = ... % ode45 limit cycle solutions
        [Xstore_wStability.stable.sol_lc2; ...
            [vec(param_i) * ones(size(Xstore_stable_lc2,1),1), Xstore_stable_lc2] ...
        ];   
    Xstore_wStability.unstable.sol_lc2 = ... 
        [Xstore_wStability.unstable.sol_lc2; ...
            [vec(param_i) * ones(size(Xstore_unstable_lc2,1),1), Xstore_unstable_lc2] ...
        ];
end