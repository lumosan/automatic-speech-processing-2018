function m = my_amdemod(s, fc, fs)
%MY_AMDEMOD Demodulate an AM signal
%   MY_AMDEMOD(S, FC, FS) Demodulates the AM signal s and returns
%   the message signal m, normalized to have values between -1 and 1. 
%   FC is the carrier frequency and FS is the sampling frequency.
%   Filters the signal with a 2nd order Butterworth filter.
%   The filtered signal has a transient at the beginning of length equal
%   to the length of the filter's impulse response, that is removed. 
order = 10;

t = 0:1/fs:(length(s)-1)/fs;
signal = s.*cos(2*pi*fc*t);

% Filter signal
[filterb, filtera] = butter(order, fc/fs);
filtered_signal = filter(filterb, filtera, signal);

% Remove transient (TODO: not sure about this)
% tfplot(filter(filterb, filtera, cat(2, 1, zeros(1, length(s) - 1))), fs,'','');
% filtered_signal = filtered_signal(20*order:length(filtered_signal));

% Remove offset
no_offset_signal = filtered_signal - mean(filtered_signal);
m = 2 * no_offset_signal;
end