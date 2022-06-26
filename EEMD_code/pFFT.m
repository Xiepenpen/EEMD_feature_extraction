function pFFT(y,Fs)
% The function has two inputs, i.e. the signal value and the sampling frequency
% The function only draws a graph, with a built-in plot, but does not use figure

%% 
t_s = 1/Fs; %sampling time steps
t_start = 0;
t_end = (length(y)-1)*t_s;
t = 0 : t_s : t_end;
Druation = t_end -t_start;  %calculate sampling time
Sampling_points = Druation/t_s +1;  %sampling points number
nfft = Sampling_points;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Frequency chart%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_f = fft(y,nfft); %Fourier transform
f_s = 1/t_s; %sampling frequency
f_x = 0:f_s/(Sampling_points -1):f_s;  %Note that this corresponds to the horizontal frequency, and the frequency resolution is f_s/(Sampling_points - 1)
plot(f_x(1:round(length(f_x)/2)),2/Sampling_points*abs(y_f(1:round(length(f_x)/2))),'k');
