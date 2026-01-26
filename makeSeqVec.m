function Xs = makeSeqVec(Xmin, Xmax, num_Xs)
% Function for producing an array of X values with inputs:
%   Xmin (lower bound), 
%   Xmax (upper bound), 
%   and num_Xs (length of vector) 
% and output:
%   Xs (the produced vector of length num_Xs).
    Xs = Xmin:((Xmax-Xmin)/num_Xs):Xmax;
end