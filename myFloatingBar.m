function bh = myFloatingBar( X, Y, BaseLine )
%myFloatingBar plots a 2-D bar graph with the zero line transposed to an
%arbitrary specified by BaseLine. Also, bars may be unevenly spaced.
%
%   Detailed explanation goes here

BarSpacing = 0;
BarColor = [0 0 0.75];

N = length(X); % Number of bars

dX = diff(X);

bh = zeros(1, N);
% Patches for 1st and Nth bar
bh(1) = patch( X(1)*[ones(1, 4)] + 0.5*dX(1)*[-1 -1 1 1] ...
    + 0.5*BarSpacing*[ 1 1 -1 -1 ], ... XData 
    BaseLine*ones(1, 4) + Y(1)*[0 1 1 0], ... YData
    BarColor);

bh(N) = patch( X(N)*[ones(1, 4)] + 0.5*dX(N-1)*[-1 -1 1 1] ...
    + 0.5*BarSpacing*[ 1 1 -1 -1 ], ... XData 
    BaseLine*ones(1, 4) + Y(N)*[0 1 1 0], ... YData
    BarColor);

for barIdx = 2:N-1
    XData = [ X(barIdx) * ones(1, 4) + ...
        [ -0.5*dX(barIdx-1)*[1 1] 0.5*dX(barIdx)*[ 1 1 ]] ...
        + 0.5*BarSpacing*[ 1 1 -1 -1 ] ];
    YData = BaseLine*ones(1, 4) + Y(barIdx)*[0 1 1 0];
        
    bh(barIdx+1) = patch( XData, YData, BarColor );
end


end

