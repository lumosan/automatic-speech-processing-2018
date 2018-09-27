function [data] = speech_signal_observation(fname, fignu, plottitle)

% This routine loads the file, plots it and returns the data
% See below the usage of this routine. It assumes that the ascii data is
% already in the file.
% data = speech_signal_observation('file name', figure number, 'plottitle')
% data is a like a variable in C.

     varname = load(fname);
     data = varname(17:length(varname));
     count = length(data);
     figure(fignu);
subplot(2, 1, 1);
     plot(data);
     title(plottitle);
     xlabel('time');
     ylabel('amplitude');
axis([0 count -30000 30000])
shift = 64;
framesize = 256;
maxframes=floor((count - framesize + shift)/shift);
energy = zeros(maxframes, 1);
for frameno=1:maxframes,
loc = (frameno-1)*shift;
sspeech = data(loc+1:loc+framesize);
energy(frameno) = sspeech' * sspeech;
end
logenergy = 10 * log10(energy);
logenergy1 = logenergy - mean(logenergy);
maxenergy = max(logenergy1);
subplot(2, 1, 2)
plot(logenergy1, 'g*');
title('Short time energy plot of the above utterance');
xlabel('time');
ylabel('Log Energy');
axis([0 maxframes 0 maxenergy])


