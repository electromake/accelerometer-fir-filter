
% Clear all local and global user-defined variables and all functions from the symbol table.
clear;

% Clear the terminal screen and move the cursor to the upper left corner.
clc;

FontSize = 20;

% Load packages 
pkg load signal
pkg load control

% g in meters per second
g = 9.81;

% Read data from CSV file
filename = "accelerometer_data_20s.txt";
data = csvread(filename);
data = data * g / 1024;

% Get X,Y,Z component of acceleration
ax = double(data(1:end,1));
ay = double(data(1:end,2));
az = double(data(1:end,3));

% Calculate magnitude of acceleration
a = sqrt(ax.^2 + ay.^2 + az.^2);

% Sampling frequency
fs = 3200;

% Nyquist frequency
fn = fs/2;

% Band pass filter frequencies [Hz]
fa =  3.0;
fb =  10.0;

% Normalized frequencies
fa_n = fa / fn; 
fb_n = fb / fn;

% Calculate filter coefficients
% Attenuation (dB)
A = 40;
% Bandwidth
BW = fb-fa;
N = (A * fs)/(22*BW)
delta1 = 1e-4;
delta2 = 1e-2;
N = 2/3 * log10(1/(delta1 * delta2)) * fs/BW

filter_order = 1829;
c = fir1(filter_order, [fa_n, fb_n], "pass");

freqz(c,1, 5000, fs)

subplot (2,1,1);
axesxy = gca;
hx = xlabel("f [Hz]");
hy = ylabel("Attenuation [dB]");
ht = title("FIR filter attenuation");
set (hx, "fontsize", FontSize);
set (hy, "fontsize", FontSize);
set (ht, "fontsize", FontSize);
axesall = findall(gcf, 'Type', 'axes');
set (axesall, "xlim", [0,40]);

subplot (2,1,2);
axesxy = gca;
hx = xlabel("f [Hz]");
hy = ylabel("Phase shift [DEG]");
ht = title("FIR filter phase shift");
set (hx, "fontsize", FontSize);
set (hy, "fontsize", FontSize);
set (ht, "fontsize", FontSize);
print -dpng -color "filter-frequency-responce.png"
clf

% Applay filter
y = filter(c, 1, az);

% Calculate time step
h = 1 / fs;

% Create time vector
t = [0:1:length(y)-1]' * h;

% Plot filtered data with grid
k = 16000;
plot(t(1:k),az(1:k),t(1:k), y(1:k), 'LineWidth',1)
hx = xlabel("t [s]");
hy = ylabel("a_z [m/s^2]");
ht = title("Accelerometer raw & filtered data");
set (hx, "fontsize", FontSize);
set (hy, "fontsize", FontSize);
set (ht, "fontsize", FontSize);
legend("Z-axes raw data", "Z-axes filtered data")
grid on;
print -dpng -color "accelerometer-Z-axes-data-2.png"
clf

% SPECTRUM OF SIGNAL
  % BEFORE FILTERING
  n = length(az);
  FT = fft(az);

  f_az = abs(FT)(1:end/2);
  arg_az = 180/pi * arg(FT)(1:end/2);
  freq_vec = [0:1:(n/2-1)] * fn/n;

  semilogx(freq_vec(2:end), 20*log10(f_az(2:end)/sqrt(n)), 'LineWidth',1)
  hx = xlabel("f [Hz]");
  hy = ylabel("Signal magnitude [dB]");
  ht = title("Accelerometer data spectrum");
  set (hx, "fontsize", FontSize);
  set (hy, "fontsize", FontSize);
  set (ht, "fontsize", FontSize);
  ylim([-100 60])
  grid on
  print -dpng -color "before-filter.png"
  clf
  
  % AFTER FILTERING
  n = length(y);
  FT = fft(y);

  f_az = abs(FT)(1:end/2);
  arg_az = 180/pi * arg(FT)(1:end/2);
  freq_vec = [0:1:(n/2-1)] * fn/n;

  semilogx(freq_vec(2:end), 20*log10(f_az(2:end)/sqrt(n)), 'LineWidth',1)
  hx = xlabel("f [Hz]");
  hy = ylabel("Signal magnitude [dB]");
  ht = title("Accelerometer data spectrum");
  set (hx, "fontsize", FontSize);
  set (hy, "fontsize", FontSize);
  set (ht, "fontsize", FontSize);
  ylim([-100 60])
  grid on
  print -dpng -color "after-filter.png"
  clf
close

% Check signal power
% sum(abs(FT).^2)/n
% sum(y.^2)