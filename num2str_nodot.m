function [ num_str ] = num2str_nodot( varargin )
%NUM2STR_NODOT Converts a number to a string but uses a delimiter other
%   than a decimal to delimit the whole number portion and the fractional
%   portion. This is especially useful for encoding decimal values in a
%   filename where the filename includes a decimal figure but it is not
%   desirable to have a decimal '.' character in the filename. The default
%   alternate delimiter is lowercase D ('d').
%
%   EXAMPLE
%
%     num2str_nodec(0.7) returns the string '0d7'.
%
%     num2str_nodec(0.7, '_') returns the string '0_7'.
%
% Created by
%   E. P. Blair
%   University of Notre Dame
%   182120RJAN12
%


num = varargin{1};

switch nargin
    case 1
        delim = 'd';
    case 2
        delim = varargin{2};
        
    otherwise
        error('NUM2STR_NODOT.M: Invalid number of input arguments.');
end % switch nargin

tempstr = num2str(num);

for char_ind = 1:length(tempstr)
    if strcmp(tempstr(char_ind), '.')
        tempstr(char_ind) = delim;
    end % END: if strcmp(tempstr(char_ind), '.')
end % END: for char_ind = 1:length(tempstr)

num_str = tempstr;

end

