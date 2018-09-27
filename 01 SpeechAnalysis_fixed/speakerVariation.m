function [] = interSpeakerDifference(file1, file2, file3, begin1, end1, begin2, end2, begin3, end3, fignu)

% This function computes the linear predictive spectrum of
% speech of different speakers and plots them in a single plot. The 
% usage of the function is as follows
%
% interSpeakerDifference('file1', 'file2', 'file3', begin1, end1, begin2, end2, begin3, end3, fignu);
%
% where file1, file2 and file3 are different speakers speech file name and
% begin1, end1, begin2, end2, begin3, end3 are the begin sample number and
% end sample number for short-time signal selection from each signal
% respectively. The sound spoken by each speaker has to be same. The linear
% prediction spectrum of all the utterances is plotted in the figure number
% fignu.

data1 = load(file1);
data2 = load(file2);
data3 = load(file3);
begin1 = begin1 + 16;
end1 = end1 + 16;
begin2 = begin2 + 16;
end2 = end2 + 16;
begin3 = begin3 + 16;
end3 = end3 + 16;
stdata1 = data1(begin1:end1);
stdata2 = data2(begin2:end2);
stdata3 = data3(begin3:end3);
lporder = 16;
lpspec1 = lpSpectrum(stdata1, lporder, 1, 512, 16000, fignu, 'lp spectrum of the file1');
%fignu = fignu + 1;
close(fignu);
lpspec2 = lpSpectrum(stdata2, lporder, 1, 512, 16000, fignu, 'lp spectrum of the file2');
%fignu = fignu +1;
close(fignu);
lpspec3 = lpSpectrum(stdata3, lporder, 1, 512, 16000, fignu, 'lp spectrum of the file3');
%fignu = fignu + 1;
close(fignu);

for n=1:512 freq(n) = ((n-1) * 16000)/(512 * 1000.0); end
figure(fignu);
plot(freq(1:256), lpspec1(1:256), 'b');
hold on
plot(freq(1:256), lpspec2(1:256), 'g');
plot(freq(1:256), lpspec3(1:256), 'r');
xlabel('Frequency (in kHz)');
ylabel('Log Amplitude');
title('Comparison of lp Spectrum, blue is for file1, green is for file2, red is for file3');
hold off


