function [] = fourierSpectrum(data, order, fignu, plottitle);

%This routine computes the fourier spectrum and plots it
%
% fourierSpectrum(data, order, fignu, 'plottitle');
%
% where data is the short time signal, order is the order of FFT
% the spectrum is plotted in a figure with figure number fignu and 
% title plottitle.

fourSpec = abs(fft(data, order));
logfourSpec = 20 * log10(fourSpec);

% assuming the sampling frequency to be 16000 the x-axis is defined
for n=1:order freq(n) = ((n-1) * 16000)/(order * 1000.0); end
figure(fignu)
     subplot(2, 1, 1)
     plot(freq, logfourSpec);
     xlabel('Frequency (in kHz)');
     ylabel('Log Amplitude');
     subplot(2, 1, 2)
     plot(freq(1:(order/2)), logfourSpec(1:(order/2)), 'r')
     title(plottitle);
     xlabel('Frequency (in kHz)');
     ylabel('Log Amplitude');
     zoom on;
%axis([0 16000/2000 20 150])


