% Matlab Code for Wireless Transmission of Data
% Dr. Ayyappadas Rajagopal

clc;
clear all;
close all;

% % Fractional number to Binary - Example usage
% fraction = -10.625; % Example negative fractional number
% disp(['Input decimal fraction: ', num2str(fraction)]);
% binaryRepresentation = fractionalToFixedPointBinary(fraction);
% disp(['Converted binary representation: ', binaryRepresentation]);
% 
% % Binary number to original fractional number - Example usage
% binaryInput = binaryRepresentation; % Example binary input
% decimalOutput = fixedPointBinaryToFraction(binaryInput);
% disp(['Reconstructed decimal fraction value: ', num2str(decimalOutput)]);




%% ---- Fractional decimal number to binary conversion ----------------- %% 
%       (16 bit Fixed point number)  
% ----------------------------------------------------------------------- %

function binaryStr = fractionalToFixedPointBinary(fraction)
    % Define fixed-point parameters
    totalBits = 16;     % Total bits
    signBits = 1;       % Sign bit
    integerBits = 9;    % Integer bits
    fractionalBits = 6; % Fractional bits

    % Determine if the number is negative
    isNegative = fraction < 0;
    fraction = abs(fraction); % Work with the absolute value

    % Separate integer and fractional parts
    integerPart = floor(fraction);
    fractionalPart = fraction - integerPart;

    % Check for overflow
    if integerPart >= 2^integerBits
        error('Integer part exceeds the maximum representable value.');
    end

    % Convert integer part to binary
    integerBinary = dec2bin(integerPart, integerBits);

    % Convert fractional part to binary
    fractionalBinary = '';
    for i = 1:fractionalBits
        fractionalPart = fractionalPart * 2;
        bit = floor(fractionalPart);
        fractionalBinary = [fractionalBinary, num2str(bit)];
        fractionalPart = fractionalPart - bit;
    end

    % Combine integer and fractional parts
    binaryStr = [integerBinary, fractionalBinary];

    % Add sign bit
    if isNegative
        % Convert to two's complement
        binaryVector = [str2num(binaryStr(:))',0]; % Append a zero for sign bit
        binaryVector = ~binaryVector; % Invert bits
        carry = 1; % Start with carry for adding 1
        for i = totalBits:-1:1
            if carry == 0
                break;
            end
            sum = binaryVector(i) + carry;
            binaryVector(i) = mod(sum, 2);
            carry = floor(sum / 2);
        end
        binaryStr = num2str(binaryVector(1:end-1)); % Exclude the sign bit
        binaryStr = strrep(binaryStr, ' ', ''); % Remove spaces
        binaryStr = ['1', binaryStr]; % Add sign bit
    else
        binaryStr = ['0', binaryStr]; % Add sign bit for positive numbers
    end
end


%% ---- Binary to Fractional decimal conversion ------------------------ %%
%       (Convert the input to binary representation)
% ----------------------------------------------------------------------- %

function decimalValue = fixedPointBinaryToFraction(binaryStr)
    % Validate input
    if length(binaryStr) ~= 16
        error('Input binary string must be 16 bits long.');
    end

    % Extract sign, integer, and fractional parts
    signBit = str2double(binaryStr(1));
    integerPartBinary = binaryStr(2:10); % Next 9 bits for integer part
    fractionalPartBinary = binaryStr(11:end); % Last 6 bits for fractional part

    % Convert binary strings to decimal
    integerPart = bin2dec(integerPartBinary);
    fractionalPart = 0;

    % Convert fractional part
    for i = 1:length(fractionalPartBinary)
        fractionalPart = fractionalPart + str2double(fractionalPartBinary(i)) * 2^(-i);
    end

    % Combine integer and fractional parts
    decimalValue = integerPart + fractionalPart;

    % Apply sign
    if signBit == 1
        % If the sign bit is 1, we need to convert to negative using two's complement
        % Calculate the maximum value for 9 bits (511) and subtract the current value
        decimalValue = decimalValue - (2^9); % Adjust for the negative range
    end
end


%% ----- QPSK Modulation ----------------------------------------------- %%
%        (Modualtion of bits to complex waveforms)
% ----------------------------------------------------------------------- %

% Parameters
fc = 1000;                 % Carrier frequency (Hz)
Fs = 10000;                % Sampling frequency (Hz)
Rb = 100;                  % Bit rate (bits/sec)
Tb = 1/Rb;                 % Bit duration (s)
bitsPerSymbol = 2;         % QPSK = 2 bits per symbol
Rs = Rb / bitsPerSymbol;   % Symbol rate (symbols/sec)
Ts = 1/Rs;                 % Symbol duration (s)

% Example binary data
data = [1  0  0  1  1  1  0  0];  % Change as needed

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
