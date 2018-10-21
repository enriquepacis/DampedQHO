function varargout = decompose(varargin)
% V = DECOMPOSE(A, B)
%     Decomposes the N-by-N matrix A into a linear combination of matrices 
%     B_i, where B_i are the elements of the cellarray B. Typically, the
%     S=(N^2)-1 generators of SU(N) are used as the set of matrices B. V is
%     a S-by-1 vector with elements defined according to Equations (2.83)
%     and (2.85) of [1]. The cellarray LAMBDA provided as the output of the
%     SU.M function is the ideal input for the B argument of this function.
%
% [V, Vo] = DECOMPOSE(A, B)
%
% EXAMPLE
%     LAMBDA = SU(4);
%     H = eye(4) - 0.025 * [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0];
%     V = decompose(H, LAMBDA);
%
% REFERENCES
%  [1] Mahler G. and V.A. Weberus, QUANTUM NETWORKS: DYNAMICS OF OPEN
%      NANOSTRUCTURES (Second Edition), Springer, Berlin, 1998. 
%
% Created by
%    Erik Blair
%    University of Notre Dame
%    010932RAUG2011
%
switch nargin
    case 1
        
    case 2
        A = varargin{1};
        B = varargin{2};
    otherwise
        error(['DECOMPOSE.M: Invalid input set. Specify an operator A', ...
            ' and a cellarray B containing a complete set of matrices', ...
            ' as a s basis for the decomposition.']);
end % END: switch nargin

szA = size(A);
n = szA(1);
s = n^2 -1;
Av = zeros(s,1);
for k = 1:s
    % B{k}
    Av(k) = trace(A*B{k});
end
A0 = trace(A);

switch nargout
    case 1
        varargout{1} = Av;
    case 2
        varargout{1} = Av;
        varargout{2} = A0;
end % END: switch nargout
        
end