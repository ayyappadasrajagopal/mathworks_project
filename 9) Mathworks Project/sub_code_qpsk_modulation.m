% Parameters
fc = 1000;                 % Carrier frequency (Hz)
Fs = 10000;                % Sampling frequency (Hz)
Rb = 100;                  % Bit rate (bits/sec)
Tb = 1/Rb;                 % Bit duration (s)
bitsPerSymbol = 2;         % QPSK = 2 bits per symbol
Rs = Rb / bitsPerSymbol;   % Symbol rate (symbols/sec)
Ts = 1/Rs;                 % Symbol duration (s)

% Example binary data
data = [1 0  0 1  1 1  0 0];  % Change as needed

% Ensure length is multiple of 2
if mod(length(data), bitsPerSymbol) ~= 0
    error('Binary data length must be multiple of 2');
end

% QPSK mapping (Gray coding)
% 00 ->  1 + j1
% 01 -> -1 + j1
% 11 -> -1 - j1
% 10 ->  1 - j1
symbols = zeros(1, length(data)/2);
for k = 1:2:length(data)
    bits = data(k:k+1);
    if isequal(bits, [0 0])
        symbols((k+1)/2) =  1 + 1j;
    elseif isequal(bits, [0 1])
        symbols((k+1)/2) = -1 + 1j;
    elseif isequal(bits, [1 1])
        symbols((k+1)/2) = -1 - 1j;
    elseif isequal(bits, [1 0])
        symbols((k+1)/2) =  1 - 1j;
    end
end

% Normalize power
symbols = symbols / sqrt(2);

% Upsample symbols to match sampling frequency
samplesPerSymbol = Fs / Rs;
I = real(upsample(symbols, samplesPerSymbol));
Q = imag(upsample(symbols, samplesPerSymbol));

% Pulse shaping (rectangular)
I = filter(ones(1, samplesPerSymbol), 1, I);
Q = filter(ones(1, samplesPerSymbol), 1, Q);

% Time vector
t = (0:length(I)-1) / Fs;

% Generate QPSK modulated carrier
qpsk_wave = I .* cos(2*pi*fc*t) - Q .* sin(2*pi*fc*t);

% Plot results
figure;
subplot(3,1,1);
stem(data, 'filled');
title('Binary Data'); xlabel('Bit index'); ylabel('Value');

subplot(3,1,2);
plot(t, I, 'b', t, Q, 'r');
title('Baseband I (blue) and Q (red) Signals'); xlabel('Time (s)');

subplot(3,1,3);
plot(t, qpsk_wave);
title('QPSK Modulated Signal'); xlabel('Time (s)'); ylabel('Amplitude');
