function [ varargout ] = getSubSysRDMIndices( varargin )
%GETSUBSYSRDMINDICES Summary of this function goes here
%   Detailed explanation goes here
%
% [RowList, ColList] = getSubSysRDMIndices( N, SubSystem )
% 
% RowList and ColList each are NScomp x NScomp x NBcomp three-dimensional
%    arrays.
%      NScomp = 2^NS
%      NBcomp = 2^NB
%
% FOR USE WITH: SubSysRDM
%
% SEE ALSO: individualRDM, getRDMindices
%

N = varargin{1};
SystemCell = sort(varargin{2});
BathCell = sort(setdiff(1:N, SystemCell));

NS = length(SystemCell);
NScomp = 2^NS;
NB = N - NS;
NBcomp = 2^NB;

RowList = zeros(NScomp, NScomp, NBcomp);
ColList = zeros(NScomp, NScomp, NBcomp);

bit = bitmatrix( N );
bitS = bitmatrix( NS );
bitB = bitmatrix( NB );
b2dS = 2.^[0:NS-1]';
b2dB = 2.^[0:NB-1]';
b2dC = 2.^[0:N-1]';
% Helper matrices MSC and MBC for converting a binary word over the system
% and the bath to a combined binary word over the closed system.
MSC = zeros(NS, N);
for syscellind = 1:NS
    MSC(syscellind, SystemCell(syscellind)) = 1;
end
MBC = zeros(NB, N);
for bathcellind = 1:NB
    MBC(bathcellind, BathCell(bathcellind)) = 1;
end

element_ct = 1;

for SSBra = 1:NScomp
    ssbrabit = bitS(SSBra, :);
        
    for SSKet = 1:NScomp
        ssketbit = bitS(SSKet, :);
        
        
        for BathStateIdx = 1:NBcomp
            bathstatebit = bitB(BathStateIdx,:);

            CompStateRow = (ssbrabit*MSC +  ...
                bathstatebit*MBC)*b2dC + 1;
            
            CompStateCol = (ssketbit*MSC +  ...
                bathstatebit*MBC)*b2dC + 1;            
            
            RowList( SSBra, SSKet, BathStateIdx ) = CompStateRow;
            ColList( SSBra, SSKet, BathStateIdx ) = CompStateCol;
        end % END: for SSBra = 1:NScomp
    end % END: for term = 1:NBcomp
end % END: for SSKet = 1:NScomp

varargout{1} = RowList;
varargout{2} = ColList;

