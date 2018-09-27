function[pdata] = select_speech(data, start, ends, fignu, plottitle);

% This routine selects the speech between the given begin and end
% points, plots it and returns it. See the usage of this routine.
%
% pdata = select_speech(data, start, end, fignu, 'plottitle');
%
% where data is the array containing the signal from which the a portion
% specified by the variable start (the starting sample number) and the
% variable end (the end sample number). The selected speech is plotted 
% in the figure fignu passed to this routine with a title passed
% through the variable plottitle. The selected region of the speech is
% returned through the variable pdata (its again like any variable in C).

pdata = data(start:ends);
figure(fignu)
     plot(pdata, 'g')
     grid on
     title(plottitle);
     xlabel('time');
     ylabel('amplitude');
     zoom on;



