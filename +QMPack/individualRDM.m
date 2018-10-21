function varargout = individualRDM(varargin) %#codegen
% INDIVIDUALRDM returns the reduced density matrix for an individual cell
% in an N-cell system. The full system density matrix must be supplied, and
% the RDM the trace is taken over the cells of non-interest.
%
% RDM = individualRDM(tgtCell, RHO)
%          Returns the reduced density matrix for the cell specified using
%          the integer tgtCell. This is the slower option, because the
%          indices over which the summation occurs in the density matrix
%          RHO must be calculated. The indices may be precalculated using
%          getRDMIndices. These, then, may simply be passed to
%          individualRDM.
%
% RDM = individualRDM(RDMIndices, RHO)
%          Returns the reduced density matrix for the cell specified using
%          the row and column indices of the full density matrix over which
%          sum is to be calculated.
%
% SEE ALSO
%          getRDMIndices.m, getSubSysRDMIndices, SubSysRDM
%
%
%
% By Enrique P. Blair
% University of Notre Dame
% Updated 181001Q NOV 2014
% 

switch nargin
    case 2
        
        % Calculate or extract the RDMIndices array
        tempArg1 = varargin{1};
        
        RHO = varargin{2};
        szRho = size(RHO);
        N = log2(szRho(1));
        
        if length(tempArg1) == 1 % fcn call: % RDM = individualRDM(tgtCell, RHO)
            tgtCell = tempArg1;
            RDMIndices = getRDMIndices(tgtCell, N);
        else % function call: RDM = individualRDM(RDMIndices, RHO)
            RDMIndices = tempArg1;
        end


    otherwise
        error('INDIVIDUALRDM.M: Invalid number of input arguments.')
end % END: switch nargin

tempRDM = zeros(2,2); % Initialize the reduced density matrix to all zeros

Nblocks = 2^(N-1);

for blockInd = 1:Nblocks
    column = RDMIndices(blockInd,:);
    c1 = column(1); c2 = column(2);
    r1 = c1; r2 = c2;
    tempRDM = tempRDM + [RHO(r1, c1) RHO(r1, c2); ...
        RHO(r2, c1) RHO(r2, c2)];
end % END: for blockInd = 1:Nblocks


switch nargout
    case 0
        varargout{1} = tempRDM;
    case 1
        varargout{1} = tempRDM;
    otherwise
        error(['INDIVIDUALRDM.M: Invalid number of output arguments. ', ...
            str2num, ' output arguments were requested.'])
end % END: switch nargout

end % END: function varargout = individualRDM(varargin)