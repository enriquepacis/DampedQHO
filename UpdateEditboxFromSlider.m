function UpdateEditboxFromSlider(varargin)
% UpdateEditboxFromSlider updates a slider linked to an EDIT box. This
% should be in the callback for the slider, and associated with the
% slider-editbox pair should be a pair of text or edit boxes which display
% the upper and lower limits of the sliding value.
%
% SYNTAX: 
%
% UpdateEditboxFromSlider( sliderHandle, textBoxHandle, ...
%         textMinValHandle, textMaxValHandle)
%
% UpdateEditboxFromSlider( sliderHandle, textBoxHandle, ...
%         textMinValHandle, textMaxValHandle, FormatString)
%
% RELATED FUNCTIONS
%   UpdateSlider, UpdateMinLim, UpdateMaxLim,
%   UpdateEditboxDiscreteFromSlider

% Required inputs
sliderhand = varargin{1};
textVal = varargin{2};
textMin = varargin{3};
textMax = varargin{4};

if nargin < 5
    FormatString = '%0.4g'; % Functional default
else
    FormatString = varargin{5};
end

xmax = str2double(get(textMax, 'String'));
xmin = str2double(get(textMin, 'String'));

slider_val = get(sliderhand, 'Value');

Dx = xmax - xmin;      dx = Dx*slider_val;
x = xmin + dx;

set(textVal, 'String', num2str(x, FormatString));
