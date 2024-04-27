clear all
close all

fs=2e6;
dec_fs=fs/150e3;

% dec_fs=fs/250000;

f = fopen('C:\Users\richa\Documents\SDR\GnuRadioprogs\MyReal', 'rb');
RealsV = fread(f, Inf, 'float');

f = fopen('C:\Users\richa\Documents\SDR\GnuRadioprogs\MyImag', 'rb');
ImagV = fread(f, Inf, 'float');

FullSig=RealsV+1i*ImagV;
FullSig = lowpass(FullSig,500e3,fs);
FT=fft(FullSig);
freq=[0:length(FT)-1]*fs/(length(FT)-1);

figure
plot(freq,fftshift(abs(FT)))

%%

phase=angle(FullSig);
Q = unwrap(phase);
m = (diff(Q));   % first derivative
%y = decimate(m,dec_fs,'fir');
% y = decimate(m,25,'fir');
y = decimate(m,16,'fir');

fsNew=fs/dec_fs;


ft_m_sig=fftshift(abs(fft(m)));
freq=[0:length(ft_m_sig)-1]*fsNew/(length(ft_m_sig)-1);
figure
plot(freq,10*log10(ft_m_sig))

fsresamp=25;
fsNew2=fsNew/fsresamp;
% y = resample(y,1,fsresamp);
desiredFs=22000;
%  [p,q] = rat(desiredFs /fsNew)
%  y = resample(y,p,q);

% [p,q] = rat(desiredFs /fsNew)
 y = resample(y,24,250);

sound((y),fsNew2)

figure
plot(m)
title('this is m')
