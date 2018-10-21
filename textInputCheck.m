function textInputCheck(textVal, textMin, textMax)
% textInputCheck should be used in the callback of an edit box with a
% numerical value stored in the 'String' property. Also associated with the
% edit box should be UI controls (either TEXT or EDIT boxes) which contain
% strings providing the upper and lower limits for the value in the edit
% box.
%
% SYNTAX:
%
% function editX_Callback(hObject, eventdata, handles)
% % hObject    handle to editPdrv (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'String') returns contents of editPdrv as text
% %        str2double(get(hObject,'String')) returns contents of editPdrv as a double
% textInputCheck(hObject, handles.textXMin, handles.textXMax);
%                  ...
%
xmin = str2double(get(textMin, 'String'));
xmax = str2double(get(textMax, 'String'));
xval = str2double(get(textVal, 'String'));

if xval > xmax
    set(textVal, 'String', num2str(xmax));
elseif xval < xmin
    set(textVal, 'String', num2str(xmin));
end
