function [] = plotSpectrogram(data, order, Window, sf, fignu, plottitle);
% this routine computes and plots the spectrogram of the given signal
% see the usage below
%
% plotSpectrogram(data, order, Window, sf, fignu, 'plottitle');
%
% where data is the array containing the signal, order is the FFT
% order, Window is the array of the window type, sf is the sampling
% frequency, fignu is the number of the plot (figure) and plottitle
% is the title of the plot.
lendata = length(data);

figure(fignu)
     subplot(2, 1, 1);
     plot(data);
     title('Speech Signal');
     xlabel('Time in sample numbers');
     ylabel('Amplitude');
     axis([1 lendata -30000 30000])
     subplot(2, 1, 2);
     specgram(data, order, sf, Window, []);
     title(plottitle);
     xlabel('Time in Seconds');
     ylabel('Frequency in Hz');


