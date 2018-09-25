function s = my_ammod(m, K, A, fc, fs)
%MY_AMMOD AM modulation of a signal
%   my_ammod(m, K, A, fc, fs) modulates the signal m using AM,
%   where fc is the carrier frequency, 
%   K and A are constants as defined in class,
%   and fs is the sampling frequency.
    t = 0:1/fs:(length(m)-1)/fs;
    s = A*(1+K*m).*cos(2*pi*fc*t);
end