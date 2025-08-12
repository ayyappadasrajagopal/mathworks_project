
% Matlab Code for Wireless Transmission of Data
% Dr. Ayyappadas Rajagopal


% Test script for fractionalToBinary function
% Define test cases
testCases = [
    struct('input', 10.625, 'precision', 5, 'expected', '1010.10100'),
    struct('input', 0.75, 'precision', 4, 'expected', '0.1100'),
    struct('input', 3.5, 'precision', 3, 'expected', '11.1'),
    struct('input', 0.1, 'precision', 10, 'expected', '0.0001100110')
];

% Run tests
for i = 1:length(testCases)
    testCase = testCases(i);
    result = fractionalToBinary(testCase.input, testCase.precision);

    % Check if the result matches the expected output
    if strcmp(result, testCase.expected)
        fprintf('Test %d passed: %s -> %s\n', i, num2str(testCase.input), result);
    else
        fprintf('Test %d failed: %s -> %s (expected: %s)\n', ...
                i, num2str(testCase.input), result, testCase.expected);
    end
end

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