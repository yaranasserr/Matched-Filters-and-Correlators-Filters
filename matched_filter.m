%1-parameters
N_bits = 1e5;        % number of bits
SNR_range = 0:2:30;  % SNR range in dB
m = input('Enter the number of samples: ');   % number of samples to represent waveform
T = input('Enter the number of sampling instant: ');   % sampling instant
amp = input('Enter the amplitude of s1: ');
s1 = amp*ones(1, m);     % rectangular signal
s2 = zeros(1, m);    % zero signal
BER = [];
CORR_BER = [];
SIM_BER = [];

%2-Generate signal
x = randi([0 1],1,N_bits);           %Vector of random values, of size 1*e5

%3-Represent each bit with proper waveform, m = 20
waveForm = zeros(1,N_bits*m);

for i = 1:length(x)   
    if x(i) == 1    %if the input bit is 1 representaion is rect of amp 1
        waveForm((i-1)*m+1:i*m) = s1;    
    else         
        waveForm((i-1)*m+1:i*m) = s2;     %else representation is zero
    end
end

 P = (1/N_bits) * sum(abs(waveForm).^2)       %This calculates the average power of the transmitted signal.
 tx_power = sum(waveForm.^2)/length(waveForm) %transmitted signal power
 
for snr = 0:2:30
    %4-Add noise to samples
    noisyWF = awgn(waveForm, snr, 'measured');
    
  %5-comparing between matched filter, correlator and simple detector:
    
  %i-Response of matched filter
    h_mf = fliplr(s1);      %reflection and shift with t=T "impulse response" of matched filter
    
    received_mf = zeros(1,N_bits);
    mf_output = zeros(1,N_bits);
    
    for i = 1 : N_bits
        noisyWF_sample = noisyWF((i-1)*m+1:i*m);    %Extracting 20 samples
        mf_sample = conv(noisyWF_sample,h_mf);
        mf_output(i) = mf_sample(T);           %convolution output at the sampling instant
    end
    
    %Calculating threshold 
    V_TH = (sum(s1.*s1)-sum(s2.*s2))/2;
    V_TH_simple = (s1+s2)/2;
    
    received_TH = zeros(1,N_bits);
    
    for j = 1:length(mf_output)
        if(mf_output(j) >= V_TH)
            received_TH(j) = 1;
        else
            received_TH(j) = 0;
        end 
    end
    
    %calculating matched filter BER
    [~, ratio] = biterr(x, received_TH);
    BER = [BER ratio]; 
     
  %ii- Correlator receiver   
     g = s1 - s2;
    corr_output = zeros(1, N_bits);
    
    for i = 1:N_bits
        start_idx = (i-1)*m+1;
        end_idx = i*m;
        corr_output(i) = sum(noisyWF(start_idx:end_idx).*g);      
    end
    
    %comparing with threshold
    corr_Received_TH = zeros(1,N_bits);
    
    for j = 1:length(corr_output)
         if(corr_output(j) >= V_TH)
            corr_Received_TH(j) = 1;
        else
            corr_Received_TH(j) = 0;
        end 
    end
    
    %calculating correlator receiver BER
    [~, ratio] = biterr(x, corr_Received_TH);
    CORR_BER = [CORR_BER ratio];
    
  %iii- simple detector
    simple_Received_TH = zeros(1,length(noisyWF));
    
    for k = 1:length(noisyWF)
         if(noisyWF(k) >= V_TH_simple)
            simple_Received_TH(k) = 1;
        else
            simple_Received_TH(k) = 0;
        end 
    end
    
    %calculating simple detector BER
    [number, ratio] = biterr(waveForm, simple_Received_TH);
    SIM_BER = [SIM_BER ratio];
end

%plotting matched filter vs correlator vs simple detector
figure;
semilogy(SNR_range, BER,'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
hold on;
semilogy(SNR_range, CORR_BER,'o','LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
hold on;
semilogy(SNR_range, SIM_BER,'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
legend('Matched filter receiver','Correlator','Simple detector');

%%
%%%%%%%%%%(7) Generalize program to have s1(t) and s2(t) as general waveforms or user defined.
%1-parameters
N_bits = 1e5;        % number of bits
SNR_range = 0:2:30;  % SNR range in dB
m = input('Enter the number of samples: ');   % number of samples to represent waveform
T = input('Enter the number of sampling instant: ');      % sampling instant
s1 = input('Enter s1: ');     
s2 = input('Enter s2: ');    
BER = [];
CORR_BER = [];
SIM_BER = [];

%2-Generate signal
x = randi([0 1],1,N_bits);           %Vector of random values, of size 1*e5

%3-Represent each bit with proper waveform, m = 20
waveForm = zeros(1,N_bits*m);

for i = 1:length(x)   
    if x(i) == 1    %if the input bit is 1 representaion is rect of amp 1
        waveForm((i-1)*m+1:i*m) = s1;    
    else         
        waveForm((i-1)*m+1:i*m) = s2;     %else representation is zero
    end
end

 P = (1/N_bits) * sum(abs(waveForm).^2)       %This calculates the average power of the transmitted signal.
 tx_power = sum(waveForm.^2)/length(waveForm) %transmitted signal power
 
for snr = 0:2:30
    %4-Add noise to samples
    noisyWF = awgn(waveForm, snr, 'measured');
    
  %5-comparing between matched filter, correlator and simple detector:
    
  %i-Response of matched filter
    h_mf = fliplr(s1);      %reflection and shift with t=T "impulse response" of matched filter
    
    received_mf = zeros(1,N_bits);
    mf_output = zeros(1,N_bits);
    
    for i = 1 : N_bits
        noisyWF_sample = noisyWF((i-1)*m+1:i*m);    %Extracting 20 samples
        mf_sample = conv(noisyWF_sample,h_mf);
        mf_output(i) = mf_sample(T);           %convolution output at the sampling instant
    end
    
    %Calculating threshold 
    V_TH = (sum(s1.*s1)-sum(s2.*s2))/2;
    V_TH_simple = (s1+s2)/2;
    
    received_TH = zeros(1,N_bits);
    
    for j = 1:length(mf_output)
        if(mf_output(j) >= V_TH)
            received_TH(j) = 1;
        else
            received_TH(j) = 0;
        end 
    end
    
    %calculating matched filter BER
    [~, ratio] = biterr(x, received_TH);
    BER = [BER ratio]; 
     
  %ii- Correlator receiver   
     g = s1 - s2;
    corr_output = zeros(1, N_bits);
    
    for i = 1:N_bits
        start_idx = (i-1)*m+1;
        end_idx = i*m;
        corr_output(i) = sum(noisyWF(start_idx:end_idx).*g);      
    end
    
    %comparing with threshold
    corr_Received_TH = zeros(1,N_bits);
    
    for j = 1:length(corr_output)
         if(corr_output(j) >= V_TH)
            corr_Received_TH(j) = 1;
        else
            corr_Received_TH(j) = 0;
        end 
    end
    
    %calculating correlator receiver BER
    [~, ratio] = biterr(x, corr_Received_TH);
    CORR_BER = [CORR_BER ratio];
    
  %iii- simple detector
    simple_Received_TH = zeros(1,length(noisyWF));
    
    for k = 1:length(noisyWF)
         if(noisyWF(k) >= V_TH_simple)
            simple_Received_TH(k) = 1;
        else
            simple_Received_TH(k) = 0;
        end 
    end
    
    %calculating simple detector BER
    [number, ratio] = biterr(waveForm, simple_Received_TH);
    SIM_BER = [SIM_BER ratio];
end

%plotting matched filter vs correlator vs simple detector
figure;
semilogy(SNR_range, BER,'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
hold on;
semilogy(SNR_range, CORR_BER,'o','LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
hold on;
semilogy(SNR_range, SIM_BER,'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
legend('Matched filter receiver','Correlator','Simple detector');