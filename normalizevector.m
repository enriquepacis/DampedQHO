function vnorm = normalizevector (V)
% NORMALIZEVECTOR normalizes a set V of N column vectors. V is passed to
%    normalizevector as a d-by-N matrix, where the N elements of V form
%    a column of the V matrix.
%
% EXAMPLES:
%
% Vnorm = normalizevector( V )
%
% E. P. Blair
% University of Notre Dame
% 241054R APR 2012
%

a = diag(sum(V.^2));
vnorm = V*(a^(-0.5));

end