function [ Lrho ] = Lindbladian( rho, L )
%Lindbladian evaluates the Lindbladian of rho given the set of Lindblad
%operators L.
%
% D = Lindbladian( rho, L )  
%
%
% E. P. Blair
% University of Notre Dame
% 311222R MAR 2014
%

[d, ~] = size(rho);

if ndims(L) == 2
    nL = 1;
else
    szL = size(L);
    nL = szL(3);
end
    
tempLrho = zeros(d);

for Lidx = 1:nL
    Lk = L(:,:,Lidx);
    Lkd = Lk';
    tempLrho = tempLrho + Lk*rho*Lkd - 0.5* anticommutator( Lkd*Lk, rho );
end % end [ for Lidx = 1:nL ]

Lrho = tempLrho;

end

