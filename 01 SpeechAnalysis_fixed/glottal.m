function []=glottal(wave,plotflag)
%---------------------------------------------------------------------------
% GLOTTAL.M
% function []=glottal(wave,plotflag)
%
% Program for estimation and plotting of diff. glottal pulse and 
% the glottal pulse. Please load the file and assign it to 'wave'.
% The differentiated glottal wave is available in the array
% 'diffglottalwave' and the glottal pulses are available in
% 'glottalwave'. Initialise 'plotflag' to enable or disable plotting.
% By default, the plots are displayed, unless the 'plotflag' is set to '0'.
%---------------------------------------------------------------------------

%------
count=size(wave,1)

% Preemphasize to remove the spetral tilt due to glottal pulse spectrum.
difwave(1) = wave(1);
difwave(2:count) = diff(wave);
if size(difwave,2)>1
   difwave = difwave';   % Convert to a column vector.
end

% Initialise the parameters.
lporder=16
nfft=512
framesize=480
shift=160
sampfreq=16000
alpha = 0.9999; % 'alpha' is the parameter for integration.
shortwin=17
thresh=40       % Threshold for residual energy.

winby2=(shortwin-1)/2;
deltafreq=sampfreq/nfft;
deltatime=1.0/sampfreq;
maxtime=framesize/sampfreq;
maxframes=floor((count-framesize+shift)/shift)
formants = zeros(maxframes,6);

% Create two arrays to plot the x-axes of the time domain and freq domain
% graphs.

for n=1:nfft/2, freq(n) = (n-1)*deltafreq/1000.0; end
for n=1:framesize, time(n) = (n-1)*deltatime*1000.0; end
win = hamming(framesize);

% Initialise the initial conditions for inverse filtering.
prevsig = zeros(lporder,1);
previntres=0;
prevint2res=0;

% Find the maximum to normalise during plotting
maxim = max(abs(wave));

echo on
% Allocating memory to speed up. Wait ...
echo off

glottalwave = zeros(count,1);
diffglottalwave = zeros(count,1);
diffspec = zeros(nfft/2,1);

echo on
% Starting frame by frame LP analysis
echo off

for frameno=1:maxframes,
   loc = (frameno-1)*shift

   % For LP analysis, store the pre-emphasised signal.
   signal = difwave(loc+1:loc+framesize);

   % For plotting, store the actual signal.
   sigbuf = wave(loc+1:loc+framesize);

   % Window the signal with a Hamming window.
   buffer = win.*signal;

   % Compute the short-time spectrum.
   %stDFT=20*log10(abs(fft(buffer,nfft)));

   % perform LP analysis.
   frameno
   arcoef = real(lpc(buffer,lporder));
   LPspect = -20*log10(abs(fft(arcoef,nfft))');   % LPspect will be a column.

   % Extract peaks of the LP spectrum.
   %diffspec(1)=0;
   %for n=2:nfft/2,diffspec(n) = LPspect(n)-LPspect(n-1); end

   diffspec(1) = 0;
   diffspec(2:nfft/2) = diff(LPspect(1:nfft/2));

   peakloc = zeros(nfft/2,1);
   nopeaks = 0;
   for n=2:nfft/2, 
      if(diffspec(n)<0 & diffspec(n-1)>=0 & n*deltafreq > 100)
         peakloc(n)=1;
         nopeaks = nopeaks+1;
      end
   end
   nopeaks

   % Store the locations of peaks.
   peakloc = find(peakloc);
   peakloc

   resocount = 1;
   if nopeaks > 3 
      newnopeaks = 3;
   else
      newnopeaks = nopeaks;
   end

   for n=1:newnopeaks,
      formants(frameno,resocount) = deltafreq * (peakloc(n)-1);
      formants(frameno,resocount+1) = LPspect(peakloc(n));
      resocount = resocount+2;
   end
      
   echo on
   % Display the formant values.
   echo off
   formants(frameno,:)

   % Compute the residual.
   augsignal = prevsig;

   augsignal(lporder+1:lporder+framesize) = signal;

   for n=lporder+1:lporder+framesize,
      predict = 0;
      for j=1:lporder,
         predict = predict + arcoef(j+1)*augsignal(n-j);
      end
      residual(n-lporder) = augsignal(n) + predict;
   end

   for n=1:lporder,
      prevsig(n) = signal(shift-lporder+n);
   end

   % Estimate the Diff glottal pulse by integrating residual of
   % pre-emphasised speech.

   intres(1)=(1-alpha)*residual(1)+alpha*previntres;
   for n=2:framesize,
      intres(n) = (1-alpha)*residual(n)+alpha*intres(n-1);
   end

   % Concatenate the diff glottal pulse.
   diffglottalwave(loc+1:loc+framesize) = intres;

   % Estimate the mean value to subtract from the diff glottal pulse.
   % This is done to remove the 'ramp' trend in the glottal wave.
   sum=mean(intres);
   shiftedintres=intres-sum;

   % Estimate the glottal pulse by integrating the dc-level shifted
   % diff glottal pulse.

   int2res(1)=shiftedintres(1)+alpha*prevint2res;
   for n=2:framesize,
      int2res(n) = shiftedintres(n)+alpha*int2res(n-1);
   end
   
   previntres = shiftedintres(shift);
   prevint2res = int2res(shift);

   % Remove the trend in the glottal pulse by subtracting regression line.
   xmean=(framesize+1)/2;
   ymean=mean(int2res);
   x2mean=(framesize+1)*(2*framesize+1)/6;

   sigmaxy=0;
   for n=1:framesize, sigmaxy=sigmaxy+n*int2res(n); end
   sigmaxy=sigmaxy/framesize;

   slope = (sigmaxy - (xmean*ymean))/(x2mean-xmean*xmean);

   for n=1:framesize, int2res(n)=int2res(n)-slope*(n-1); end

   % Concatenate the glottal pulse.
   for n=1:framesize,
      glottalwave(n+loc) = int2res(n);
   end

   % Compute short time energy of residual to identify peaks
   for n=1:winby2, tempres(n)=0; end
   for n=winby2+1:winby2+framesize, tempres(n)=residual(n-winby2); end
   for n=winby2+framesize+1:framesize+2*winby2, tempres(n)=0; end

   for n=1:framesize, 
      resenergy(n)=0;
      index=n+winby2;
      for k=-winby2:winby2,
         resenergy(n)=resenergy(n)+(tempres(index+k)*tempres(index+k));
      end
   end
   
   maxenergy=max(resenergy);
   % Median filter the energy
   %medfilt1(resenergy,5);

   % Threshold the energy in each frame as 'thresh' precent 
   % of the peak value.
   if resenergy(1) > ((thresh/100)*maxenergy)
      resenergy(1) = 1;
   else
      resenergy(1) = 0;
   end

   peakcount=0;
   for n=2:framesize, 
      if resenergy(n) > ((thresh/100)*maxenergy)
         resenergy(n) = 1;
      else
         resenergy(n) = 0;
      end

      if (resenergy(n-1)==0) & (resenergy(n)==1)
         peakcount=peakcount+1;
      end
   end
   peakcount

   % Find the negative maxima in the diff. glottal pulse which is
   % not dc-level shifted.
   peakcounter=0;
   m=2;
   flag=0;
   while (peakcounter<=peakcount) & (m<=framesize),
      if(resenergy(m-1)==0) & (resenergy(m)==1) 
         flag=1;

      else m=m+1;
      end
      k=1;
      if flag==1
         while flag==1 & m<framesize,
            peakregions(k)=intres(m);
            k=k+1;
            m=m+1;
            if(resenergy(m)==0) flag=0; end
         end
         peakcounter=peakcounter+1;
         peaks(peakcounter)=min(peakregions);
         %if peakcounter==peakcount break, end
      end   
      %if peakcounter==peakcount break, end
   end
 
   echo on
   % Display the peak values of diff. glottal pulse
   echo off
   peaks(1:peakcount)

   % Display the plots, if 'plotflag' is not 0.
   if plotflag~=0
      % Plot the actual signal
      subplot(321),plot(time,sigbuf,'y')
      grid
      xlabel('time (ms)')
      ylabel('pressure')
      axis([0 maxtime*1000 -maxim maxim])

      % Plot the short-time spectrum
      %subplot(322),plot(freq,stDFT(1:nfft/2),'g')
      %grid
      %xlabel('freq (kHz)')
      %ylabel('dB')
      %axis([0 sampfreq/2000 40 230])

      % Plot the LP spectrum
      subplot(322),plot(freq,LPspect(1:nfft/2),'c')
      grid
      xlabel('freq (kHz)')
      ylabel('dB')
      axis([0 sampfreq/2000 -70 100])

      % Plot the glottal pulse.
      subplot(323),plot(time,int2res,'c')
      grid
      xlabel('time (ms)')
      ylabel('pulse')
      axis([0 maxtime*1000 -0.5 0.5])

      % Plot the residual of the pre-emphasised speech.
      subplot(324),plot(time,residual,'g')
      grid
      xlabel('time (ms)')
      ylabel('residual')
      axis([0 maxtime*1000 -maxim/2.5 maxim/2.5])

      % Plot the differentiated glottal pulse.
      subplot(325),plot(time,intres,'y')
      grid
      xlabel('time (ms)')
      ylabel('diff pulse')
      axis([0 maxtime*1000 -0.08 0.08])

      % Plot the thresholded residual short time energy of the 
      % pre-emphasised speech.
      subplot(326),plot(time,resenergy,'c')
      grid
      xlabel('time (ms)')
      ylabel('thresholded energy')
      axis([0 maxtime*1000 0 1.1])

      pause(2)
   end
end

% Plot all the results obtained.
figure(2)

% Plot the speech signal.
subplot(311),plot(wave,'c')
grid
xlabel('time')
ylabel('pressure')

% Plot glottal pulses
subplot(312),plot(glottalwave+ones(length(glottalwave),1)*8,'g')
grid
xlabel('time')
ylabel('Volume vel')

% Plot diff glottal pulses
subplot(313),plot(diffglottalwave,'y')
grid
xlabel('time')
ylabel('diff Volume vel')

figure(3)

% Plot 1st formant contour
subplot(311),plot(formants(:,1),'g')
grid
xlabel('time')
ylabel('First formant')

% Plot 2nd formant contour
subplot(312),plot(formants(:,3),'c')
grid
xlabel('time')
ylabel('Second formant')

% Plot 3rd formant contour
subplot(313),plot(formants(:,5),'y')
grid
xlabel('time')
ylabel('Third formant')

figure(3)

% Plot the first formant energy.
subplot(616),plot(formants(:,2),'c')
grid
xlabel('time')
ylabel('F1 energy (dB)')

% Plot the second formant energy.
subplot(615),plot(formants(:,4),'w')
grid
xlabel('time')
ylabel('F2 energy (dB)')

% Plot the third formant energy.
subplot(614),plot(formants(:,6),'b')
grid
xlabel('time')
ylabel('F3 energy (dB)')

% Plot 1st formant contour
subplot(613),plot(formants(:,1),'g')
grid
xlabel('time')
ylabel('First formant')

% Plot 2nd formant contour
subplot(612),plot(formants(:,3),'c')
grid
xlabel('time')
ylabel('Second formant')

% Plot 3rd formant contour
subplot(611),plot(formants(:,5),'y')
grid
xlabel('time')
ylabel('Third formant')

     figure(4)
     specgram(wave, 256, 16000, hamming(256), []);
     xlabel('Time in Seconds');
     ylabel('Frequency in Hz');
     title('Spectrogram of the signal');
