function UpdateSlider(sliderhand, textVal, textMin, textMax)
% UpdateSlider updates a slider which is tied to an EDIT box.
%
% Use UpdateSlider in the callback for an EDIT box when the slider and edit
% box should be linked. There also should be TEXT/EDIT boxes specifying the
% upper limit of the editable variable.
%
% Also, textInputCheck should occur before calling UpdateSlider.
%
% Recommended syntax:
%    
%    function myEditBox_X_CallBack(hObject, eventdata, handles)
%    % hObject    handle to editPdrv (see GCBO)
%    % eventdata  reserved - to be defined in a future version of MATLAB
%    % handles    structure with handles and user data (see GUIDATA)
%
%    textInputCheck(hObject, handles.textXMin, handles.textXMax);  
%
%    UpdateSlider(handles.slider_X, hObject, handles.textXMin, handles.textXMax)
%
xmin = str2double(get(textMin, 'String'));
xmax = str2double(get(textMax, 'String'));
xval = str2double(get(textVal, 'String'));

Dx = xmax - xmin; dx = xval - xmin;

set(sliderhand, 'Value', dx/Dx);
