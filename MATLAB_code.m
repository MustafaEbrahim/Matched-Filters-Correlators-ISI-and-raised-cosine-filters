close all;
clear;
clc
%% control flags
pulse_width = 1; % Ts =  1 s
pulse_samples =  5 ; % No. of samples to represent the shaping pulse signal
imp_tr_length = 10; %no of the impulses in the sampling
Spacing = pulse_width/pulse_samples; % To de-normalize the time index y[n]= y[nTs]
%% generate the transmitted signal
pulse = [5 4 3 2 1]/sqrt(55); %generate the unit energy pulse-shaping signal
imp_tr = randi([0 1],1,imp_tr_length);  % generate the random logic bit-stream
imp_tr_mapped = 2*imp_tr-1; % map logic 1 into (+1) & logic 0 into (-1)
imp_tr_upsampled = upsample(imp_tr_mapped,pulse_samples); %create the sampling period by upsampling (Ts= 1s ,Tp = 20 ms --> 5 samples)
imp_tr_upsampled = imp_tr_upsampled(1:end-4) ;%to elimiate the spacing after the last pulse
tx_out = conv(pulse,imp_tr_upsampled); % produce the tranmitted signal
%% plot the pulse signal and its PAM
figure();
subplot(2,2,1);
pulse_x = (0:length(pulse)-1).*Spacing;
plot(pulse_x,pulse);
xlim([0 1])
xlabel("Time (Sec)");
ylabel("Amplitude")
hold on;
stem(pulse_x,pulse);
title("Pulse Shaping Function");
subplot(2,2,2);
imp_x = (0:length(imp_tr_upsampled)-1).*Spacing;
stem(imp_x,imp_tr_upsampled);
xlabel("Time (Sec)");
ylabel("Amplitude")
title("PAM Sampling Function");
subplot(2,2,[3 4]);
tx_out_x = (0:length(tx_out)-1).*Spacing;
plot(tx_out_x,tx_out);
hold on;
stem(imp_x,imp_tr_upsampled);
set(gca,"XAxisLocation",'origin');
title("Transmitted Signal");
sgtitle("Pulse Shaping Scheme @ Tx");
figure;
plot(tx_out_x,tx_out);
set(gca,'XAxisLocation','origin');
xlabel("Time (Sec)");
ylabel("Amplitude")
title("Output of the TX");
%% Matched Filter 
tr_x = (1:imp_tr_length)-Spacing;
MF_filter_TF = fliplr(pulse);
MF_filter_out = conv(tx_out,MF_filter_TF);
figure()
subplot(2,1,1);
MF_x = (0:length(MF_filter_out)-1).*Spacing;
plot(MF_x,MF_filter_out);
hold on;
stem(tr_x,imp_tr_mapped);
set(gca,"XAxisLocation",'origin');
title("Rx signal using Matched filter");
xlabel("Time (Sec)");
ylabel("Amplitude")
%% Non-matched Filter
NMF_filter_TF = ones(1,pulse_samples)/sqrt(length(pulse));
NMF_filter_out = conv(tx_out,NMF_filter_TF);
subplot(2,1,2);
plot(MF_x,NMF_filter_out,'k');
hold on;
stem(tr_x,imp_tr_mapped);
set(gca,"XAxisLocation",'origin');
title("Rx signal using Non-Matched filter");
xlabel("Time (Sec)");
ylabel("Amplitude")
%% Matching filter using correlator
corr_out = correlate(tx_out,pulse); %p(t)*p(t)
% To show the correlator & the matched filter on the same plot
figure()
subplot(2,1,1);
stem(tr_x,imp_tr_mapped);
title("Transmitted PCM sequence");
subplot(2,1,2);
plot(MF_x,MF_filter_out,'r');
hold on;
plot(MF_x(1:end-4),corr_out,'b');
hold on;
stem(tr_x,corr_out(5:5:end),'g');
title("Correlator output vs Matched filter output")
xlabel("Time (Sec)");
ylabel("Amplitude")
set(gca,'XAxisLocation','origin');
xlim([0 tr_x(end)]);
legend("Matched filter output","Correlator output ","Correlator output sampled",'Location','best');
%% Noise Analysis
% Generate the binary sequency with numerous bits to introduce the Pe as BER
imp_tr_length = 10000;
imp_tr = randi([0 1],1,imp_tr_length);  % generate the random logic bit-stream
imp_tr_mapped = 2*imp_tr-1; % map logic 1 into (+1) & logic 0 into (-1)
imp_tr_upsampled = upsample(imp_tr_mapped,pulse_samples); %create the sampling period by upsampling (Ts= 1s ,Tp = 20 ms --> 5 samples)
imp_tr_upsampled = imp_tr_upsampled(1:end-4) ;%to elimiate the spacing after the last pulse
tx_out = conv(pulse,imp_tr_upsampled); % produce the tranmitted signal
% Generate the noise 
AWGN = randn(1,length(tx_out)); % ~N(0,1) Scaled the variance by scaling the standard deviation 
% create 2 arrays to get the practical BER for each case
k=1;
BER_MF_practical = zeros(1,8);
BER_NMF_practical = zeros(1,8);
BER_theoritical =  zeros(1,8);
Eb = 1;
for i= logspace(-0.2,0.5,8)
%Since we are using unit-energy p(t). Hence, 1/N0 E [-2,5]dB
    N0 = Eb/i;
    noise = AWGN * sqrt(N0/2);
% A new signal due to the addition of the AWGN Noise @ the reciever "Z(t)" 
z = noise + tx_out;
% Passing the signal (z(t)) to the matched/non-matched filter                            
MF_filter_at_Rx = correlate(z,pulse);
NMF_filter_at_Rx = correlate(z,NMF_filter_TF);
% Sampling the output @ the reciever
MF_filter_at_Rx = receiver_sampler(MF_filter_at_Rx,pulse_samples);
NMF_filter_at_Rx = receiver_sampler(NMF_filter_at_Rx,pulse_samples);
                        % Reciever Decision 
%Sample the filter output then take a threshold @ 0 :
% If >0 --> Hence, the tranmitted symbol is p(t) --> 1
% If <0 --> Hence , the transmitted symbol is -p(t) --> 0
MF_Rx_decision = decision_block(MF_filter_at_Rx,0);
NMF_Rx_decision = decision_block(NMF_filter_at_Rx,0);
% Calculate the BER by just counting the mismatches 
% between the actual PCM & decision upon the constellation
 % We can compare directly with the mapping bec it is deterministic
check_MF = (MF_Rx_decision ~= imp_tr_mapped);
check_NMF = (NMF_Rx_decision ~= imp_tr_mapped);
BER_MF_practical(k) = length(check_MF(check_MF==1))/imp_tr_length;
BER_NMF_practical(k) = length(check_NMF(check_NMF==1))/imp_tr_length;
BER_theoritical(k) = 0.5*erfc(sqrt(i));
k=k+1;
end
% plot the BER for Matched/Unmatched/Theoritical filters on semi-log axis
figure();
BER_x = -2:5;
semilogy(BER_x,BER_theoritical,'g');
hold on;
grid  on ;
semilogy(BER_x,BER_MF_practical,'b');
hold on;
semilogy(BER_x,BER_NMF_practical,'r');
legend("Theoritical","Matched Filter",...
        "Non-Matched Filter");
xlabel("Eb/N0");
ylabel("BER");
title("BER across Matching / Non-matching Filter");
%% ISI & Raised Cosine
imp_tr_length = 100;
imp_tr = randi([0 1],1,imp_tr_length);  % generate the random logic bit-stream
imp_tr_mapped = 2*imp_tr-1; % map logic 1 into (+1) & logic 0 into (-1)
imp_tr_upsampled = upsample(imp_tr_mapped,pulse_samples);
roll_off = [0 1];
delay = [2 8];
for i=1:2
    for j=1:2
       srrc_tx = rcosine(pulse_width,pulse_samples,'sqrt',roll_off(i),delay(j));
       %        srrc_tx = rcosdesign(roll_off(i),delay(j),pulse_samples,'sqrt');
       tx_out = conv(imp_tr_upsampled,srrc_tx,'valid'); %To exclude the last zeros from the convolution
       % Assume the channel is white & noise-free
       % Suppose using matching srrc filter
       srrc_Rx = fliplr(srrc_tx);
       figure
       impz(srrc_tx); 
       title(['Roll-off = ',num2str(roll_off(i)),', Delay = ',num2str(delay(j))]);
       % Knowing that the Raised Cosine versions are symmetric --> Srrc_Rx = Srrc_Tx
       Rx_in = conv(tx_out,srrc_Rx,'valid');
       eyediagram(tx_out,2*pulse_samples);
       title(['At Tx , Roll-off = ',num2str(roll_off(i)),', Delay = ',num2str(delay(j))]);
       eyediagram(Rx_in(1:end-6),2*pulse_samples);
       title(['At Rx , Roll-off = ',num2str(roll_off(i)),', Delay = ',num2str(delay(j))]);
    end
end
function Z_signal = correlate (Y_signal,X_signal)
Z_signal = zeros(1,length(Y_signal));
for i=1:length(X_signal):length(Y_signal)
corr_window = zeros(1,length(X_signal));
    for j= 0:length(corr_window)-1
    corr_window(j+1) = Y_signal(i+j);
    Z_signal(i+j) = sum(corr_window.*X_signal);
    end 
end
end
function z = receiver_sampler (x,n) % x ---> output of the filter , n ---> Sampling period
z = zeros(1,floor(length(x)/n));
for i = 1:length(z)
   z(i) = x(n*i); 
end
end
function Y_Signal = decision_block(X_Signal , threshold) 
    Y_Signal = zeros (1, length (X_Signal));
    Y_Signal(X_Signal > threshold) = 1;
    Y_Signal(X_Signal <= threshold) = -1;
end 