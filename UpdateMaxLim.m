function UpdateMaxLim(textVal, textMin, textMax, sliderhand)
%  UpdateMaxLim should be used in the callback of an edit box specifying
%  the upper limit of a variable in a four-object UI control system
%  consisting of an edit box for a variable X coupled to a slider.
%  There also should be edit boxes for the minimum and maximum values of
%  the variable X.
%
%  SYNTAX
%
%  >> UpdateMaxLim(textVal, textMin, textMax, sliderhand)
%
%  DEPENDENCIES
%     UpdateSlider
%     
%  SEE ALSO:
%     UpdateMinLim, 
%
%  E.P. Blair
%  University of Notre Dame
%  Updated: 230124Q NOV 2014
xmin = str2double(get(textMin, 'String'));
xmax = str2double(get(textMax, 'String'));
xval = str2double(get(textVal, 'String'));

% xval is assumed to be greater than or equal to xmin
if xmax < xval
    set(textMax, 'String', num2str(xval));
end

UpdateSlider(sliderhand, textVal, textMin, textMax);
