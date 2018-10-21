function varargout = purity ( varargin )
% PURITY returns the purity of a density matrix RHO.
%
% P = purity(RHO)
%

rho = varargin{1};

P = trace(rho*rho);

varargout{1} = P;

end