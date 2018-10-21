function varargout = entropyShannon( varargin )
% ENTROPYSHANNON calculates the Shannon entropy of a density matrix
%      RHO or a probability distribution P (a row or a column vector).
%
% SSHN = entropyShannon ( RHO ) calculates the Shannon entropy
%      SSHN for RHO.
%
% By E. P. Blair
%    University of Notre Dame
%    161314R FEB 2012
%

rho = varargin{1};

[a, b] = size(rho);
if (a == b) && a > 1 % RHO IS A DENSITY MATRIX
    q = real(diag(rho));
    p = q(q~=0);
elseif ismember(1, [a, b]) % RHO IS A COLUMN (VECTOR)
    p = rho(rho~=0);
end

Sshn = -sum(p.*log2(p));

varargout{1} = Sshn;