function vocalize( speechstring )
% This function uses the voice synthesizer to speak the string SPEECHSTRING
% on a Mac or a PC.
%
% E. Blair
% 215212R APR 2012
% University of Notre Dame
% 
if ismac
    system(['say -v Alex ', speechstring])
elseif ispc
    tts(speechstring)
end

end
