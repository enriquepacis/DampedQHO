function A = recompose(A0, Ak, Lambda)
% A = RECOMPOSE(A0, Ak, LAMBDA)
%  Recomposes matrix A from linearly combining the basis matrices Lambda
%     in a weighted sum with coefficients A0 and Ak. This implements
%     Equation (2.83) in [1]. V is a vector of length s=n^2-1; Lambda is a
%     complete set of s matrices; A0 is a scalar.
%
%  This was designed for reconstructing a density matrix rho from a
%  coherence vector lambda consisting of lambda_k = lambda_1 ... lambda_n.
%  Thus, A0 maps to
%      lambda_0 = trace(rho) = 1
%  and Ak maps to the coherence vector.
%  LAMBDA maps to the generators of SU(n).
%  Thus, use the syntax:
%    > rho = recompose(1, lambda, LAMBDA)
%  where 1 is a constant, lambda is the coherence vector, and LAMBDA are
%  the s = 2n^2-1 generators SU(n).
%
% REFERENCES
%  [1] Mahler G. and V.A. Weberus, QUANTUM NETWORKS: DYNAMICS OF OPEN
%      NANOSTRUCTURES (Second Edition), Springer, Berlin, 1998. 
%
% Created by
%    Erik Blair
%    University of Notre Dame
%    010958RAUG2011
%

s = length(Ak);
n = sqrt(s+1);

A = (1/n)*A0*eye(n);
for k = 1:s
    A = A + 0.5*Ak(k)*Lambda{k};
end
 
end % END: function A = recompose(A0, Ak, Lambda)