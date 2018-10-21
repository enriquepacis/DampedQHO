function [ rho_r_X, rho_r_Y ] = BipartiteTraceout( rho, dX, dY )
% BipartiteTraceout returns the reduced density matrix of the two
% components in a bipartite system XY
%
%   [ rho_r_X, rho_r_Y ] = BipartiteTraceout( rho, dX, dY )
%
%   INPUTS
%   ======
%
%   rho      Density matrix for the global system XY.
%
%   dX       The dimension of subsystem X.
%
%   dY       The dimension of subsystem Y.
%
%   Assume that the global system is initially in a product state
%
%         rho(0) = kron( rhoX(0), rhoY(0) ) ,
%
%   and that it undergoes unitary evolution under the quantum Liouville
%   equation:
%
%         rho(t) = U(H,t) * ( rho(0) ) * (U(H,t)') ,
%
%   where U(H, t) = expm( - ( 1i/hbar ) * H * t).
%
%

rho_r_X = zeros(dX);
rho_r_Y = zeros(dY);

% Trace rho over the degrees of freedom of Y to get rho_r_X
for k = 1:dY
    Xaddend = rho(k:dY:k+(dX-1)*dY, k:dY:k+(dX-1)*dY);
    rho_r_X = rho_r_X + Xaddend;
end % for k = 1:dY

% Trace rho over the degrees of freedom of X to get rho_r_Y
for k = 1:dX
    Yaddend = rho(dY*k-(dY-1):dY*k, dY*k-(dY-1):dY*k );
    rho_r_Y = rho_r_Y + Yaddend;
end % for k = 1:dY


end

