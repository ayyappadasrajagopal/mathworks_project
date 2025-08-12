% Matlab Code for Wireless Transmission of Data
% Dr. Ayyappadas Rajagopal

% Example usage
% fraction = 10.625; % Example fractional number
% precision = 5;     % Number of bits for the fractional part
% binaryRepresentation = fractionalToBinary(fraction, precision);
% disp(['Binary representation: ', binaryRepresentation]);
% 
% 
% % Define state-space matrices
% A = [0 1; -5 -2];  % State matrix
% B = [0; 3];        % Input matrix
% C = [0 1];        % Output matrix
% D = [0];          % Feedthrough matrix
% 
% % Create the state-space model
% sys = ss(A, B, C, D, 'InputDelay', 0);  % Discrete-time model
% 
% % Define simulation parameters
% t = 0:0.1:10;  % Time vector
% u = sin(t);    % Input signal (e.g., a sine wave)
% 
% % Simulate the system response
% [y, t] = lsim(sys, u, t);
% 
% % Display the output measurements
% disp('Output Measurements:');
% disp(y);
% 
% % Convert measurements to binary
% precision = 10;  % Number of bits for fractional part
% binaryMeasurements = arrayfun(@(x) fractionalToBinary(x, precision), y, 'UniformOutput', false);
% 
% % Display binary measurements
% disp('Binary Measurements:');
% disp(binaryMeasurements);



% Example usage
fraction = -10.625; % Example negative fractional number
precision = 5;      % Number of bits for the fractional part
binaryRepresentation = fractionalToBinary(fraction, precision);
disp(['Binary representation: ', binaryRepresentation]);

function binaryStr = fractionalToBinary(fraction, precision)
    % Validate input
    if ~isnumeric(fraction)
        error('Input must be a numeric value.');
    end

    % Initialize variables
    isNegative = fraction < 0;  % Check if the number is negative
    fraction = abs(fraction);    % Work with the absolute value for conversion
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

    % If the original number was negative, convert to two's complement
    if isNegative
        % Convert binary string to a number
        totalBits = length(binaryStr) - 1; % Exclude the dot
        % Create a binary vector
        binaryVector = zeros(1, totalBits);
        for i = 1:totalBits
            binaryVector(i) = str2double(binaryStr(i));
        end
        % Invert bits
        binaryVector = ~binaryVector;
        % Add 1 to the least significant bit
        carry = 1;
        for i = totalBits:-1:1
            if carry == 0
                break;
            end
            sum = binaryVector(i) + carry;
            binaryVector(i) = mod(sum, 2);
            carry = floor(sum / 2);
        end
        % Convert back to string
        binaryStr = num2str(binaryVector);
        binaryStr = strrep(binaryStr, ' ', ''); % Remove spaces
    end
end



% Example usage
fraction = -10.625; % Example negative fractional number
precision = 5;      % Number of bits for the fractional part
binaryRepresentation = fractionalToBinary(fraction, precision);
disp(['Binary representation: ', binaryRepresentation]);
