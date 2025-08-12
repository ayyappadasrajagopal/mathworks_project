% Matlab Code for Wireless Transmission of Data
% Dr. Ayyappadas Rajagopal

clc;
clear all;
close all;

% Fractional number to Binary - Example usage
fraction = -10.625; % Example negative fractional number
disp(['Input decimal fraction: ', num2str(fraction)]);
binaryRepresentation = fractionalToFixedPointBinary(fraction);
disp(['Converted binary representation: ', binaryRepresentation]);

% Binary number to original fractional number - Example usage
binaryInput = binaryRepresentation; % Example binary input
decimalOutput = fixedPointBinaryToFraction(binaryInput);
disp(['Reconstructed decimal fraction value: ', num2str(decimalOutput)]);

% Example usage
binaryInput = [0 0 0 1 1 1 1 0]'; % Example binary input
carrierFrequency = 1000; % Carrier frequency in Hz
sampleRate = 8000; % Sample rate in Hz
qpskModulation(binaryInput, carrierFrequency, sampleRate);


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
% (Here the modualtion of bits to corresponding complex waveforms is done)
% ----------------------------------------------------------------------- %

function qpskModulation(binaryData, carrierFreq, sampleRate)
    % Ensure binaryData is a row vector
    if iscolumn(binaryData)
        binaryData = binaryData';
    end

    % Step 1: Convert binary data to symbols (2 bits per symbol)
    numSymbols = length(binaryData) / 2;
    symbols = zeros(1, numSymbols);

    for i = 1:numSymbols
        bits = binaryData(2*i-1:2*i);
        % Map bits to QPSK symbols
        if isequal(bits, [0 0])
            symbols(i) = 1; % Phase 0
        elseif isequal(bits, [0 1])
            symbols(i) = 1i; % Phase 90 degrees
        elseif isequal(bits, [1 1])
            symbols(i) = -1; % Phase 180 degrees
        elseif isequal(bits, [1 0])
            symbols(i) = -1i; % Phase 270 degrees
        end
    end

    % Step 2: Create time vector for the carrier wave
    t = (0:1/sampleRate:(numSymbols-1)/sampleRate);

    % Step 3: Modulate the symbols onto the carrier wave
    carrierWave = cos(2 * pi * carrierFreq * t);
    modulatedSignal = real(symbols * exp(1i * 2 * pi * carrierFreq * t));

    % Plot the results
    figure;
    subplot(3, 1, 1);
    stem(binaryData, 'filled');
    title('Input Binary Data');
    xlabel('Bit Index');
    ylabel('Bit Value');

    subplot(3, 1, 2);
    plot(real(symbols), imag(symbols), 'o');
    title('QPSK Constellation Diagram');
    xlabel('In-Phase');
    ylabel('Quadrature');

    subplot(3, 1, 3);
    plot(t, modulatedSignal);
    title('Modulated Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

