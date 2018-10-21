function [ drhodt ] = Lindblad( rho, hbar, H, L )
%Lindblad evalautes the instantantaneous time-derivative of the density
%    matrix rho using the Lindblad equation. This is valid for constant H.
%
%   [ drhodt ] = Lindblad( rho, hbar, H, L )
%
% E. P. Blair
% University of Notre Dame
% 311939R MAR 2014
%

drhodt = -(1i/hbar)*commutator(H, rho) + Lindbladian(rho, L);

end

