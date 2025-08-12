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

function binaryStr = fractionalToFixedPointBinary(fraction)
    % Define fixed-point parameters
    totalBits = 16;  % Total bits
    signBits = 1;    % Sign bit
    integerBits = 9; % Integer bits
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





