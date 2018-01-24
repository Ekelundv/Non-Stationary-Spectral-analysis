%% Sinusoid signal

clear all;

% Generating data.

fs = 1000;   % Sampling frequency
t = (0:fs)/fs;  % Time.
%A = [1 2 3 4 5];  % Amplitudes.
A = [1 2];  % Amplitudes.
%f = [0.1; 0.4; 0.08; 0.2; 0.07]; % % Sinusoid frequencies (column vector)
f = [0.1; 0.01;]; % % Sinusoid frequencies (column vector)
d = length(f); % Model order.
xn = A(1)*sin(2*pi*f(1)*t) + A(2)*sin(2*pi*f(2)*t) + 10*randn(size(t));
m = 20; % Size of the Covariance.

% Producing the estimated periodogram

[Pxx,F] = periodogram(xn,hamming(length(xn)),length(xn),fs);
plot(F,10*log10(Pxx))
xlabel('Hz')
ylabel('dB')
title('Modified Periodogram Power Spectral Density Estimate')

order = sinorder(xn,1,fs);

% Estimate the frequencies. Order the estimates.
fESPRIT = sort( esprit(xn, d, m)/2/pi )'
fWSF    = sort( MODE( covM(xn,m), d ) )'
fRELAX  = sort( relax(xn,d)/2/pi )
%fMUSIC  = sort( rootmusic(xn,d)/2/pi )';

%f = [0.1; 0.4; 0.08; 0.2; 0.07];

%% rational model AR(8) - process
clear all;
N= 800;
A = [1 -0.98 0.89 -0.89 0.85 -0.97 0.89 -0.9 0.91];
C = [1];
e = rand(N,1);
y = filter(C,A,e);
plot(y)
%%
figure;
subplot(2,1,1)
plot(y)
title('Simulated data using an AR(8)-process with 800 data points')
subplot(2,1,2)
periodogram(y,[],4096)

order = sinorder(y,var(y),N);
%% MODEL ESTIMATION
N = 800;
k=1;
AR1=arx(y,k);
yp = predict(AR1,y)
BIC = N*log(var(yp)) + 1*log(N)

%%
N = 800;
k=2;
AR2=arx(y,k);
yp = predict(AR2,y)
BIC = N*log(var(yp)) + k*log(N)
%%
N = 800;
k=3;
AR3=arx(y,k);
yp = predict(AR3,y)
BIC = N*log(var(yp)) + 1*log(N)
%%
N = 800;
k=4;
AR4=arx(y,k);
yp = predict(AR4,y)
BIC = N*log(var(yp)) + 1*log(N)
%%
N = 800;
k=5;
AR5=arx(y,k);
yp = predict(AR5,y)
BIC = N*log(var(yp)) + 1*log(N)
%%
N = 800;
k=6;
AR6=arx(y,k);
yp = predict(AR6,y)
BIC = N*log(var(yp)) + 1*log(N)
%%
N = 800;
k=7;
AR7=arx(y,k);
yp = predict(AR7,y)
BIC = N*log(var(yp)) + 1*log(N)
%%
N = 800;
k=8;
AR8=arx(y,k);
yp = predict(AR8,y)
BIC = N*log(var(yp)) + 1*log(N)
%% MATLAB EXAMPLE - MODIFIED
clear all;
fs = 1000;                % Sampling frequency
t = (0:fs)/fs;            % One second worth of samples
A = [1 2 3 4 5];                % Sinusoid amplitudes (row vector)
f = [140;90; 150; 160; 80];            % Sinusoid frequencies (column vector)
e = 1*randn(size(t));
xn = A*sin(2*pi*f*t) + e;
% The three last lines are equivalent to
% xn = sin(2*pi*150*t) + 2*sin(2*pi*140*t) + 0.1*randn(size(t));
SNR = var(A*sin(2*pi*f*t))/var(e)

% The periodogram estimate of the PSD can be computed using |periodogram|.
% In this case, the data vector is multiplied by a Hamming window to
% produce a modified periodogram.

[Pxx,F] = periodogram(xn,[],length(xn),fs);
plot(10*log10(Pxx))
xlabel('Hz')
ylabel('dB')
title('Modified Periodogram Power Spectral Density Estimate')

%% A. Jokobsson example
clear all;
N   = 64;        % Number of samples.
f   = [.09 .1];  % Frequencies
amp = [1 2];     % Amplitudes
m   = 20;        % Size of the covariance matrix. 
P   = 1024;      % Zeropadding used by the periodogram.

d = length(f);
x = amp(1)*exp( i*f(1)*2*pi*(1:N)' + pi*rand*i ) + amp(2)*exp( i*f(2)*2*pi*(1:N)' + pi*rand*i );
w = .1*( randn(N,1) + i*randn(N,1) )/sqrt(2);
y = x + w;

% Estimate the frequencies. Order the estimates.
fMUSIC  = sort( rootmusic(y,d)/2/pi )';
fESPRIT = sort( esprit(y, d, m)/2/pi )';
fWSF    = sort( MODE( covM(y,m), d ) )';
fRELAX  = sort( relax(y,d)/2/pi );

% Estimate the periodogram
Y = fftshift( abs( fft( y, P )/N ) ).^2;
ff = (0:P-1)/P-.5;
semilogy( ff, Y );
fPER = findpeaks(Y,d); % Finds the d most dominat peaks
fPER = sort( ff(fPER) );

disp( sprintf('True frequencies:  [ %f, %f ]', f(1), f(2)) )
disp( sprintf('  MUSIC            [ %f, %f ], rMSE %f ', fMUSIC(1), fMUSIC(2), sqrt(sum(abs(f-fMUSIC).^2)) ) )
disp( sprintf('  ESPRIT           [ %f, %f ], rMSE %f', fESPRIT(1), fESPRIT(2), sqrt(sum(abs(f-fESPRIT).^2))) )
disp( sprintf('  WSF              [ %f, %f ], rMSE %f', fWSF(1), fWSF(2), sqrt(sum(abs(f-fWSF).^2))) )
disp( sprintf('  RELAX            [ %f, %f ], rMSE %f', fRELAX(1), fRELAX(2), sqrt(sum(abs(f-fRELAX).^2))) )
disp( sprintf('  Periodogram      [ %f, %f ], rMSE %f\n', fPER(1), fPER(2), sqrt(sum(abs(f-fPER).^2))) )
disp( sprintf('Periodogram resolution: 1/N = %f (roughly)', 1/N ) )   % On average, the periodogram has a resolution of about 1/N.
disp( sprintf('Grid resolution limit:  1/P = %f \n', 1/P ) )          % On average, a grid search cannot beat this resolution.

%%
clear all;
fs = 1000; % Number of samples.
t = (1:fs)/fs;
f = [0.09; 0.1]; % Frequencies. 
amp = [1 2]; % Amplitudes
e = 1*randn(1,fs); % Noise, increase/decrease VAR with 1

x = amp*sin(2*pi*f*t) + e;    
[Pxx,F] = periodogram(x,[],length(x),fs);
plot(10*log10(Pxx))

%% 
clear all;
% Generating some data

N=1024;
fs=200;
f = [0.09 0.1]; % Frekvenser
d = length(f);
A = [1 2];
ts= 1/fs;
t = ts*(0:N-1); 
y =  A(1)*sin(2*pi*f(1)*t) + A(2)*sin(2*pi*f(2)*t) + 0.1*randn(size(t));
plot(t,y)
m = 20;

% Estimating the peaks
fESPRIT = sort( esprit(y, d, m)/2/pi )'
fWSF    = sort( MODE( covM(y,m), d ) )'
fRELAX  = sort( relax(y,d)/2/pi )
