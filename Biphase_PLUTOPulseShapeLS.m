clear all;
close all;




% Message to be transmitted
msg = 'hello My name is Richard';
binary = dec2bin(msg);
str = bin2dec(binary);

sps =8; % Number of samples per symbol (oversampling factor)

% Convert message to binary
msg_bin = dec2bin(double(msg), 8);
save_msg_bin =msg_bin;

[row cols]=size(msg_bin);
msg_symbolsBiphse = reshape(msg_bin,[row*cols 1]);

% Reshape binary message into 4-bit symbols
msg_symbols = reshape(msg_bin.', 4, []).';
save_msg_symbols =str2num(msg_symbolsBiphse);


% Convert binary symbols to decimal values
msg_decimal = bin2dec(msg_symbols);

filtlen = 1; % Filter length in symbols
rolloff = 0.75; % Filter rolloff factor
%
rrcFilter = rcosdesign(rolloff,filtlen,sps);

msg_symbolsBiphse=str2num(msg_symbolsBiphse);

for ii=1:length(msg_symbolsBiphse)
    if(msg_symbolsBiphse(ii)<1)
        msg_symbolsBiphse(ii)=-1;
    end
end

txFiltSignal = upfirdn(msg_symbolsBiphse,rrcFilter,sps,1);
% figure(2)
% plot(1:length(txFiltSignal),real(txFiltSignal),'-b',1:length(msg_symbolsBiphse),msg_symbolsBiphse,'--*r')
% title('Filtered Tx sig')



figure(50)
plot(real(txFiltSignal),'-.')






%% =========================================================
%===========================================================


fs =1e6;
dt=1/fs;
Pw=500e-6;
Bw=0.1e6;
K=Bw/Pw;
f1=1e4;
chiplen=1e-5;
Fc=1500.0e6;
SampsChip=round(chiplen*fs);
%Pcode=[1 1 1 1 1 0 0 1 1 0 1 0 1];

%
% t=[0:dt:chiplen*length(Pcode)-dt];
%
% PcodeSamp=[];
% for ii=1:length(Pcode)
%   PcodeSamp=[PcodeSamp Pcode(ii)*ones(1,SampsChip)];
% end
%
% figure(1)
% plot(PcodeSamp)
% txWaveform=exp(1i*pi*(f1*t+PcodeSamp));
t=txFiltSignal;
%txWaveform=hilbert(txSig);
txWaveform=hilbert(txFiltSignal);
txWav=hilbert(txFiltSignal);
% txWaveform=hilbert(PcodeSamp);

figure(2)
plot(1:length(txWaveform), real(txWaveform),'b',1:length(txWaveform),imag(txWaveform),'--r')

% figure(30)
% plot(real(txWaveform),imag(txWaveform),'*')
% xlabel('real inphase')
% ylabel('Imag quadrature')




% %%add zeros to make it pulse
Lzero=10;
%  txWaveform=[zeros(1,Lzero*length(t))+0i*zeros(1,Lzero*length(t)) transpose(txSig) zeros(1,Lzero*length(t))+0i*zeros(1,Lzero*length(t)) ];
%
txWaveform=[zeros(1,Lzero*length(t))+0i*zeros(1,Lzero*length(t)) transpose(txWav) zeros(1,Lzero*length(t))+0i*zeros(1,Lzero*length(t)) 0 ];
%

figure(3)
plot(1:length(txWaveform),real(txWaveform),'-b',1:length(txWaveform),real(txWaveform),'r.')
title('Waveform to be sent to PLUTO')
txWaveform=transpose(txWaveform);

%
% MF=xcorr(txWaveform);
% figure(4)
% plot(abs(MF))
% title('Waveform to be sent to PLUTO Autocorrelation')

% scatterplot(txWaveform)
pause(1)
tx = sdrtx('Pluto');
tx.CenterFrequency = Fc;
tx.BasebandSampleRate = fs;
tx.Gain = -20;
transmitRepeat(tx,txWaveform);



rx = sdrrx('Pluto',...
    'BasebandSampleRate',fs,'CenterFrequency',Fc,'SamplesPerFrame',1*length(txWaveform));
%rx = sdrrx('Pluto')
framesToCollect=8;
frameSize=rx.SamplesPerFrame;
%% Template 1
% Perform data collection then offline processing
data = zeros(frameSize, framesToCollect);
% Collect all frames in continuity
for frame = 1:framesToCollect
    [d,valid,of] = rx();

    % Collect data without overflow and is valid
    %     if ~valid
    %         warning('Data invalid')
    %     elseif of
    %         warning('Overflow occurred')
    %     else
    data(:,frame) = d;
    %     end
end
%
% % Process new live data
% sa1 = dsp.SpectrumAnalyzer;
% for frame = 1:framesToCollect
%     sa1(data(:,frame)); % Algorithm processing
% end
% pause(1)


figure(10)
imagesc(real(data))

figure(9)
plot((real(data(:,1))))
title('example frame')


rxFiltSignal = upfirdn(data(:,1),rrcFilter,1,sps); % Downsample and filter
rxFiltSignal = rxFiltSignal(filtlen + 1:end - filtlen); % Account for delay




figure(19)
plot(real(rxFiltSignal ))
title('example frame')


pause(3)
Myrelease;

[match1]=conv(txWav,data(:,1));

sampcoord=round(1*[1:length(match1)]);
figure(101)
plot(abs(match1))


%%

SigInfo=real(rxFiltSignal);

info_indices=find(abs(SigInfo)>1000);

min_info=min(info_indices);
max_info=max(info_indices);

figure;
plot(SigInfo(min_info:max_info))

NewSig=rxFiltSignal(min_info:max_info);

ampest_I=norm(real(NewSig))/norm(real(msg_symbolsBiphse));

ampest_Q=norm(imag(NewSig))/norm(imag((msg_symbolsBiphse)));


%% Scatterplots

% scatplot = scatterplot(sqrt(sps)*...
%     data(:,1),...
%     sps,0);
% hold on;
% % scatterplot(rxFiltSignal(1:5e3),1,0,'bx',scatplot);
% scatterplot(rxFiltSignal,1,0,'bx',scatplot);
% title('Received Signal, Before and After Filtering');
% legend('Before Filtering','After Filtering');
% % axis([-5 5 -5 5]); % Set axis ranges
% hold off;



%% Demodulate QAM

SigInfoScaled_I=real(NewSig)/ampest_I;
SigInfoScaled_Q=imag(NewSig)/ampest_Q;



% overI=find(SigInfoScaled_I>0);
% SigInfoScaled_I(overI)=3;
% overI=find(SigInfoScaled_I<-3);
% SigInfoScaled_I(overI)=-3;
%
% overQ=find(SigInfoScaled_Q>3);
% SigInfoScaled_Q(overQ)=3;
% overQ=find(SigInfoScaled_Q<-3);
% SigInfoScaled_Q(overQ)=-3;


FullSig=SigInfoScaled_I+1i*SigInfoScaled_Q;

figure(5);
plot(1:length(FullSig),round(real(FullSig)),'-*')
RxPluto=real(FullSig);

RxPluto=round(real(FullSig(2:end)));
RxPluto2=round(real(FullSig));
bit0=find(RxPluto<0);
RxPluto(bit0)=0;
bit0=find(RxPluto2<0);
RxPluto2(bit0)=0;

RxPlutoShift=circshift(RxPluto2,length(RxPluto2)-1);
figure(22);plot(1:length(save_msg_symbols),save_msg_symbols,'*-r',1:length(RxPlutoShift),RxPlutoShift,'-ob')

figure(23);plot(1:length(save_msg_symbols),save_msg_symbols,'*-r',1:length(RxPluto2),RxPluto2,'-ob')

figure(6)
plot(1:length(save_msg_symbols),  save_msg_symbols,'*-r',1:length(RxPluto),RxPluto,'--ob')

if(mod(length(RxPluto),8)==0)
    RxPlutoUse=RxPluto;
end
if(mod(length(RxPluto2),8)==0)
    RxPlutoUse=RxPluto2;

    %     RxPlutoUse=RxPlutoShift(1:end-1);
end

Rxsigclean=dec2bin((  RxPlutoUse));


Rx_bits=reshape(Rxsigclean,[],8);

demodulated_msg = char(bin2dec(Rx_bits));

% Display original and demodulated messages
disp(['Original message: ', msg]);
% disp(['Demodulated message: ', demodulated_msg]);
disp(['Demodulated message: ',transpose(demodulated_msg)]);







%  demodulated = qamdemod(round(FullSig), 16);
% %
% % % Convert decimal values to binary symbols
%
% demodulated_symbols = dec2bin(demodulated, 4);
% % % Reshape binary symbols into binary message
%  demodulated_bin = reshape(demodulated_symbols.', [], 1);
% offset=mod(length(demodulated_bin),8);
%
% %
% % % Convert binary message to ASCII characters
% demodulated_msg = char(bin2dec(reshape(demodulated_bin(offset+1:end), 8, []).'));
%
%
%
% % % Display original and demodulated messages
%  disp(['Original message: ', msg]);
%
% disp(['Demodulated message: ',transpose(demodulated_msg)]);
%
%
%
