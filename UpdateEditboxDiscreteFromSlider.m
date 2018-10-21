function UpdateEditboxDiscreteFromSlider(varargin)
% UpdateEditboxFromSlider( sliderHandle, textBoxHandle, ...
%         textMinValHandle, textMaxValHandle)
%
% UpdateEditBoxFromSlider( sliderHandle, textBoxHandle, ...
%         textMinValHandle, textMaxValHandle, FormatString)
%
% SUGGESTED SYNTAX:
%
%   % --- Executes on slider movement.
%   function sliderMySlider_Callback(hObject, eventdata, handles)
%   % hObject    handle to sliderTimeIndex (see GCBO)
%   % eventdata  reserved - to be defined in a future version of MATLAB
%   % handles    structure with handles and user data (see GUIDATA)
%   
%   % Hints: get(hObject,'Value') returns position of slider
%   %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%
%   UpdateEditboxFromSlider( hObject, handles.editValue, ...
%           handles.textValueMin, handles.textValueMax )
%

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

Dx = xmax - xmin;
x = xmin + round(slider_val*(Dx));

set(textVal, 'String', num2str(x, FormatString));
