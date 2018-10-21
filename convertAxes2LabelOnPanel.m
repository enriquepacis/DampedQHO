function convertAxes2LabelOnPanel( varargin )
% convertAxes2Label( AxHand, AxString, ... ) makes an AXES object look
%      like the background upon which it is set, and it writes the string
%      AxString in the axes.
%
%      Options: the 'FontName' and 'FontSize' properties for the text may be
%      set.

AxesHandle = varargin{1}; % Required input
AxesString = varargin{2};

% defaults
Interpreter = 'latex';
FontSz = 12;
FontName = 'Helvetica';

% FUNCTIONAL ARGUMENTS
args = varargin(3:end);
while length(args) >= 2
    prop = args{1};
    val = args{2};
    args = args(3:end);
    switch prop
        case 'FontSize'
            FontSz = val;
        case 'FontName'
            FontName = val;
        case 'Interpreter'
            Interpreter = val;
    end % END [ switch prop ]
end % END [ while length(args) >= 2 ]

PanelHandle = get(AxesHandle, 'Parent');
PanelColor = get(PanelHandle, 'BackgroundColor');
%Blk = [0 0 0];
cla(AxesHandle);
set(AxesHandle, 'Color', PanelColor, ...
    'XTickLabel', [], 'XTick', [], 'XColor', PanelColor, ...
    'YTickLabel', [], 'YTick', [], 'YColor', PanelColor );

axes( AxesHandle );
myText = text(0, 0, AxesString);
set(myText, 'Interpreter', Interpreter, ...
    'FontSize', FontSz, 'FontName', FontName);