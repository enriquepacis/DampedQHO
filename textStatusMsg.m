function textStatusMsg( targettextbox, textString, MsgColor)
% textStatusMsg (MyTextBox, TextString, TextColor)
%

set(targettextbox, 'String', [datestr(now), ': ', textString], ...
    'ForegroundColor', MsgColor);

drawnow;

