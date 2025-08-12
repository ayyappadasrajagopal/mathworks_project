% Matlab Code for Wireless Transmission of Data
% Dr. Ayyappadas Rajagopal

function binaryStr = fractionalToBinary(fraction, precision)
    % Validate input
    if fraction < 0
        error('Input must be a non-negative fractional number.');
    end

    % Initialize variables
    integerPart = floor(fraction);
    fractionalPart = fraction - integerPart;
    binaryStr = '';

    % Convert integer part to binary
    binaryStr = [binaryStr, dec2bin(integerPart), '.'];

    % Convert fractional part to binary
    while precision > 0
        fractionalPart = fractionalPart * 2;
        bit = floor(fractionalPart);
        binaryStr = [binaryStr, num2str(bit)];
        fractionalPart = fractionalPart - bit;
        precision = precision - 1;
    end
end

% Example usage
fraction = 10.625; % Example fractional number
precision = 5;     % Number of bits for the fractional part
binaryRepresentation = fractionalToBinary(fraction, precision);
disp(['Binary representation: ', binaryRepresentation]);