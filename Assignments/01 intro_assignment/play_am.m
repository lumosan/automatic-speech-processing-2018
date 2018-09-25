function play_am()
%PLAY_AM Loads and demodulates a signal
%   The signal is loaded and demodulated, then is downsampled by some
%   factor and is reproduced.
    file_name = 'am_signal.mat';
    fs = 100000; % Hz
    fc = 20000; % Hz
    downsampling_factor = 10;
    
    file = load(file_name);
    modulated_m = file.am_signal;
    
    demodulated_m = sol_amdemod(modulated_m, fc, fs);

    downsampled_m = downsample(demodulated_m, downsampling_factor);
    downsampled_fs = fs / downsampling_factor;
     
    sound(downsampled_m, downsampled_fs);
end