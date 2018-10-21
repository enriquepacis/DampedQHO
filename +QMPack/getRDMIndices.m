function varargout = getRDMIndices( varargin )
%GETRDMINDICES Returns indices for the RDM for the kTH cell of an N-cell
% array.
%
% For an N-cell system of two-dot QCA cells, the output of GETRDMINDICES
% will be a 2^(N-1)-by-2 array, where the first column contains row
% indices and the second column contains column indices.
%
% The output of this function will be used by the function INDIVIDUALRDM.
%
% EXAMPLE:
%
% RDMIndices = getRDMIndices(TargetCell, N)
%
%
% FOR USE WITH: individualRDM
%
% SEE ALSO: getSubSysRDMIndices, SybSysRDM
%
%

switch nargin
    case 2
        tgtCell = varargin{1};
        N = varargin{2};
    otherwise
        error(['GETRDMINDICES.M: Invalid number of inputs (', ...
            num2str(nargin), ' inputs were provided).']);
end % END: switch nargin

Nblocks = 2^(N-1);
% disp(['Target Cell is: ', num2str(tgtCell)])
RDMIndices = zeros(Nblocks, 2);
for BlockInd = 1:Nblocks
    tempStr = fliplr(dec2bin(BlockInd-1, N-1));
    c1str = [tempStr(1:tgtCell-1), '0', ...
            tempStr(tgtCell:end)];
    c1 = bin2dec(fliplr(c1str)) + 1;
    c2str = [tempStr(1:tgtCell-1), '1', ...
            tempStr(tgtCell:end)];
    c2 = bin2dec(fliplr(c2str)) + 1;
    RDMIndices(BlockInd, :) = [c1 c2];
end % END: for BlockInd = 1:Nblocks

switch nargout
    case 0
        varargout{1} = RDMIndices;
    case 1
        varargout{1} = RDMIndices;
end % END: switch nargout

end

