function varargout = outerprod(v1, v2);
% P = outerprod(v1, v2)
%   Matrix P is the projection |v1> <v2| (or outer product) of two n-by-1
%          kets represented by column vectors |v1> and |v2>.
%
% Created by Erik Blair
% University of Notre Dame
% July 19, 2011
%

P = v1*v2';

switch nargout
    case 0
        varargout{1} = P;
    case 1        
        varargout{1} = P;
    otherwise
        error('OUTERPROD.M: too many input arguments.');
end