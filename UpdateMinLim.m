function UpdateMinLim(textVal, textMin, textMax, sliderhand)
%  UpdateMinLim should be used in the callback of an edit box specifying
%  the upper limit of a variable in a four-object UI control system
%  consisting of an edit box for a variable X coupled to a slider.
%  There also should be edit boxes for the minimum and maximum values of
%  the variable X.
%
%  SYNTAX
%
%
%  DEPENDENCIES
%     UpdateSlider
%     
%  SEE ALSO:
%     UpdateMaxnLim, xmin = str2double(get(textMin, 'String'));
%
xmax = str2double(get(textMax, 'String'));
xval = str2double(get(textVal, 'String'));

% xval is assumed to be greater than or equal to xmin
if xmin > xval
    set(textMin, 'String', num2str(xval));
end

UpdateSlider(sliderhand, textVal, textMin, textMax);
