function [ drho_v_dt ] = LindbladN( t, rho0, H, L, hbar )
%LindbladN applies the Lindblad equation with a Hamiltonian H and a set of 
%   NL Lindblad operators to evolve a density matrix rho0.
%
%  SYNTAX
%
%      [ t, rho_t ] = ode45(@LindbladN, ...
%              [0 tmax], rho_v_0, [], H, L, hbar );
% 
%   rho_v_0   Initial density matrix for an N-dimensional system, provided
%             as a (N^2)-by-1 vector input argument by means of
%             transforming the square density matrix rho0 into rho_v_0 via
%
%             >>  rho_v_0 = rho0(:);
%
%   H         The constant Hamiltonian   
%
%   L         NL N-by-N Lindblad operators for the system, packaged in a
%             [N N NL] three-dimensional MATLAB array
%
%   rho_t     nt-by-N^2 time-varying density matrix. Each row corresponds
%             to the density matrix at time sample point.
%
% E. P. Blair
% University of Notre Dame
% 212048R JAN 2014
%

N = sqrt(length(rho0));
rho = reshape(rho0, [N N]);

drhodt = Lindblad( rho, hbar, H, L );

drho_v_dt = drhodt(:);
