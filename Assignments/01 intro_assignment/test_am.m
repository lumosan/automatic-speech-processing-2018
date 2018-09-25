function test_am()
%TEST_AM Create and modulate a test signal
%   The signal is created and modulated.
%   Original and modulated signals are plot.
    clear all;
    close all;
    
    f_info = 10; % Hz
    fs = 4000; % Hz
    fc = 300; % Hz
    
    K = 1; A = K;
    d = 1; % signal duration in s
    
    t = linspace(0, d, d*fs + 1);
    m = 1/2 * cos(2*pi*f_info*t);
    
    % Plot 1: original signal
    tfplot(m, fs, 'm_{am}', 'Message signal')
    
    % Modulate
    modulated_m = my_ammod(m, K, A, fc, fs);
    
    % Plot 2: modulated signal
    tfplot(modulated_m, fs, 's_{am}', 'AM modulated signal')

    % Demodulate
    demodulated_m = my_amdemod(modulated_m, fc, fs);
    
    % Plot 3: demodulated signal
    tfplot(demodulated_m, fs, 'm_{am}(est)', 'Recovered message signal')
    
    % Reproduce sounds
    sound(m, fs)
    sound(demodulated_m, fs)
end

