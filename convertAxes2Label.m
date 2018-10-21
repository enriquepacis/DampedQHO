function convertAxes2Label( varargin )
% convertAxes2Label( AxHand, AxString, ... ) makes an AXES object look
%      like the background upon which it is set, and it writes the string
%      AxString in the axes.
%
%      Options:
%
%          'FontName'
%          'FontSize'
%          'Rotation'
%          'HorizontalAlignment'
%          'VerticalAlignment'
%          'Rotation'
%
% This function is best used

AxesHandle = varargin{1}; % Required input
AxesString = varargin{2};

% defaults
Interpreter = 'latex';
FontSz = 12;
FontName = 'Helvetica';
Rotation = 0;
VerticalAlignment = 'middle';
HorizontalAlignment = 'center';

% FUNCTIONAL ARGUMENTS
args = varargin(3:end);
while length(args) >= 2
    prop = args{1};
    val = args{2};
    args = args(3:end);
    switch prop
        case 'FontSize'
            if ischar(val)
                FontSz = str2double(val);
            else
                FontSz = val;
            end
        case 'FontName'
            FontName = val;
        case 'Interpreter'
            Interpreter = val;
            
        case 'HorizontalAlignment'
            HorizontalAlignment = val;

        case 'VerticalAlignment'
            VerticalAlignment = val;

        case 'Rotation'
            if ischar(val)
                Rotation = str2double(val);
            else
                Rotation = val;
            end
            TextRotation = Rotation
    end % END [ switch prop ]
end % END [ while length(args) >= 2 ]


FigColor = get(gcf, 'Color');
%Blk = [0 0 0];
cla(AxesHandle);
set(AxesHandle, 'Color', FigColor, ...
    'XTickLabel', [], 'XTick', [], 'XColor', FigColor, ...
    'YTickLabel', [], 'YTick', [], 'YColor', FigColor );

axes( AxesHandle );
myText = text(0, 0, AxesString);
set(myText, 'Interpreter', Interpreter, ...
    'FontSize', FontSz, 'FontName', FontName, ...
    'Rotation', Rotation, ...
    'HorizontalAlignment', HorizontalAlignment, ...
    'VerticalAlignment', VerticalAlignment);
