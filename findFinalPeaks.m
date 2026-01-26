function [minX, maxX, isLimit, Xstore] = findFinalPeaks(y, var_num, bi_frac, tolLim, Xstore, param_i, ic)
% Function for finding the lower and upper bounds of the limiting cycle, if
% it exists, and for storing the solution appropriately with inputs:
%   y (the time series solutions), 
%   var_num (integer corresponding to the state variable under investigation: 1 for N, 2 for P, 3 for Z),
%   bi_frac (the fraction of time-points considered as burn-in), 
%   tolLim (the threshold within which limit cycles are deemed to be steady),
%   Xstore (the long-term solutions for compartment X given an array of bifurcation parameter values),
%   param_i (the index for the bifurcation parameter value),
%   and ic (index of initial condition).
% and outputs:
%   minX (the minimum value of the limiting cyle), 
%   maxX (the maximum value of the limiting cycle), 
%   isLimit (1 if the limits appear to correspond to a limiting cycle, 2 if not),
%   and Xstore (the updated long-term solutions for compartment X).

    % Find limits
    num_tps = size(y,1);                                        % Number of time-points
    index_after_bi = ceil(bi_frac .* num_tps);
    maximaX = sort(findpeaks(y(index_after_bi:end,var_num)));   % All local peaks in the solution post burn-in
    minimaX = sort(-findpeaks(-y(index_after_bi:end,var_num))); % All local troughs in the solution post burn-in      

    if ~isempty(maximaX) && ~isempty(minimaX) ...           % MAXIMA AND MINIMA EXIST
        && ~isscalar(maximaX) && ~isscalar(minimaX)
        minX = min(minimaX); 
        maxX = max(maximaX);

        if abs(maxX-minX)<tolLim                            % Oscillations appear to be a stationary point (difference b/n max and min very small)
            isLimit = 1;                                    
        elseif abs(maximaX(end-1)-maxX)<tolLim && ...       % Oscillations appear to be a limit cycle
                abs(minimaX(2)-minX)<tolLim                 
            isLimit = 2;                                    
        else                                                % Oscillations do not appear to reach limitting pattern
            isLimit = 0;                                    
        end
    else                                                    % MAXIMA AND MINIMA DO NOT EXIST
        minX = min(y(index_after_bi:end,var_num));
        maxX = max(y(index_after_bi:end,var_num));

        % if abs(maxX-minX)<tolLim                            % Monotonous approach to stationary point (difference b/n max and min very small)
        %     isLimit = 1;                                      
        % else                                                % No stationary point or limit is approached
        %     isLimit = 0;                                    
        % end    
        isLimit = 0;
    end

    % Store ODE45 solution and type of limiting behaviour
    if isLimit == 1                     % Approximately stationary solution
        Xstore(param_i,ic,3) = mean([minX,maxX]);
        Xstore(param_i,ic,5) = Inf;     % Indicating true Limit, with NaN col 4 indicating no limit cycle
    elseif isLimit == 2                 % Limit cycle
        Xstore(param_i,ic,3) = minX; 
        Xstore(param_i,ic,4) = maxX;
        Xstore(param_i,ic,5) = Inf;     % Indicating true Limit
    elseif isLimit == 0                 % Limit not reached
        Xstore(param_i,ic,3) = minX; 
        Xstore(param_i,ic,4) = maxX;
        Xstore(param_i,ic,5) = NaN;     % Indicating unknown true limit
    end

    % Presume stability for ODE45 solution
    if isLimit == 0                     % If no ODE45 solution
        Xstore(param_i,ic,6) = NaN;
    elseif Xstore(param_i,ic,2) == 0    % If the FSOLVE solution is unstable
        Xstore(param_i,ic,6) = 1;       % ...the ODE45 solution is stable
    elseif Xstore(param_i,ic,2) == 1    % If the FSOLVE solution is stable
        Xstore(param_i,ic,6) = 0;       % ...the ODE45 solution is unstable
    end
end