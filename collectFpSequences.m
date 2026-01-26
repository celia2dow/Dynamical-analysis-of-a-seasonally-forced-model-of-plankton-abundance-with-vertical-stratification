function [Xruns, nRuns] = collectFpSequences(Xstore_wStability_field, tolX)
% Function for collecting sequences of stationary fixed points in a run
% with inputs:
%   Xstore_wStability_field (specific field containing ordered values of long-term solutions given a particular stability),
%   and tolX (the threshold difference beyond which it is deemed that a jump has occured).
% and outputs:
%   Xruns (a cell array containing the collected sequences of stationary fixed points),
%   and nRuns (the number of runs found in this field).

    if ~isempty(Xstore_wStability_field)
        sortedField = sortrows(Xstore_wStability_field);
        endsWithoutJumps = [false; find(abs(diff(sortedField(:,2)))>tolX)] + 1;
        endsWithoutJumps = [endsWithoutJumps; size(Xstore_wStability_field,1) + 1];
        nRuns = length(endsWithoutJumps) - 1;
        Xruns = cell(nRuns,1);
        index_lo = endsWithoutJumps(1:end-1);
        index_up = endsWithoutJumps(2:end) - 1;
        for k = 1:nRuns
            Xruns{k} = sortedField(index_lo(k):index_up(k),:);
        end
    else
        Xruns = [];
        nRuns = 0;
    end
end