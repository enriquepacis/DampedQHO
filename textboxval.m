function val = textboxval( varargin )
%textboxval combines STRING2DOUBLE with GET. This allows the user to code
%
%        x = textboxval( hgObject )
%
%    instead of
%
%        x = str2double( get(hgObject, 'String') );
%
%  SEE ALSO
%     settextval.m
%
%  E. P. Blair
%  031418R APR 2014
%  University of Notre Dame
%

hgObject = varargin{1};

val = str2double( get(hgObject, 'String') );

end

