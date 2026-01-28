function [Xstore_stable, Xstore_unstable] = findStabUnstabSols(Xstore, tol)
% Function for sorting a given type of solution and extracting the unique
% stable and unstable long-term solutions and their stabilities with 
% inputs:
%   Xstore (the long-term solutions for compartment X given an array of bifurcation parameter values),
%   and tol ( tolerance for the solution being steady in the long-term).
% and outputs:
%   Xstore_stable (array storing unique, ordered values of stable long-term solutions),
%   and Xstore_unstable (array storing unique, ordered values of stable long-term solutions).
    
    % Sort steady state solutions for given bifurcation parameter in terms of ODE45 solution value
    Xsorted = sortrows(Xstore);    
    num_sol_cols = size(Xsorted,2)-1;

    % Find IC indices of stable and unstable solutions
    stableXics = squeeze(Xsorted(:,end)) == 1;   
    unstableXics = squeeze(Xsorted(:,end)) == 0;

    % Find values of stable and unstable solutions
    stableXvals = squeeze(Xsorted(stableXics,1:end-1));   
    unstableXvals = squeeze(Xsorted(unstableXics,1:end-1));  
   
    % Number of distinct solutions for this bifurcation parameter value
    if any(stableXics) % Number of stable groups
        % Indices of beginning and end of consecutive sequence of nearby solutions
        endsXstableGrps = [false; find(abs(diff(stableXvals(:,1)))>tol)] + 1;  
        endsXstableGrps = [endsXstableGrps; length(stableXvals(:,1)) + 1];

        nStabGrps = length(endsXstableGrps) - 1;
        icStab_l = endsXstableGrps(1:end-1);
        icStab_u = endsXstableGrps(2:end) - 1;
        Xstore_stable =  NaN(nStabGrps, num_sol_cols); 
        for k = 1:nStabGrps
            Xstore_stable(k,:) = ...
                mean(stableXvals(icStab_l(k):icStab_u(k),:));
        end
    else
        Xstore_stable = [];
    end

    if any(unstableXics) % Number of unstable groups   
        % Indices of beginning and end of consecutive sequence of nearby solutions
        endsXunstableGrps = [false; find(abs(diff(unstableXvals(:,1)))>tol)] +1; 
        endsXunstableGrps = [endsXunstableGrps; length(unstableXvals(:,1)) + 1];
   
        nUnstabGrps = length(endsXunstableGrps) - 1;  
        icUnstab_l = endsXunstableGrps(1:end-1);
        icUnstab_u = endsXunstableGrps(2:end) - 1;
        Xstore_unstable =  NaN(nUnstabGrps, num_sol_cols); 
        for k = 1:nUnstabGrps 
            Xstore_unstable(k,:) = ...
                mean(unstableXvals(icUnstab_l(k):icUnstab_u(k),:));
        end
    else
        Xstore_unstable = [];
    end
end