% Matlab Code for Wireless Transmission of Data
% Dr. Ayyappadas Rajagopal


% function binaryStr = fractionalToBinary(fraction, precision)
%     % Validate input
%     if ~isnumeric(fraction)
%         error('Input must be a numeric value.');
%     end
% 
%     % Check for NaN
%     if isnan(fraction)
%         binaryStr = 'NaN'; % Return 'NaN' as a string for NaN input
%         return;
%     end
% 
%     % Initialize variables
%     isNegative = fraction < 0;  % Check if the number is negative
%     fraction = abs(fraction);    % Work with the absolute value for conversion
%     integerPart = floor(fraction);
%     fractionalPart = fraction - integerPart;
%     binaryStr = '';
% 
%     % Convert integer part to binary
%     binaryStr = [binaryStr, dec2bin(integerPart), '.'];
% 
%     % Convert fractional part to binary
%     while precision > 0
%         fractionalPart = fractionalPart * 2;
%         bit = floor(fractionalPart);
%         binaryStr = [binaryStr, num2str(bit)];
%         fractionalPart = fractionalPart - bit;
%         precision = precision - 1;
%     end
% 
%     % If the original number was negative, convert to two's complement
%     if isNegative
%         % Convert binary string to a number
%         totalBits = length(binaryStr) - 1; % Exclude the dot
%         % Create a binary vector
%         binaryVector = zeros(1, totalBits);
%         for i = 1:totalBits
%             binaryVector(i) = str2double(binaryStr(i));
%         end
%         % Invert bits
%         binaryVector = ~binaryVector;
%         % Add 1 to the least significant bit
%         carry = 1;
%         for i = totalBits:-1:1
%             if carry == 0
%                 break;
%             end
%             sum = binaryVector(i) + carry;
%             binaryVector(i) = mod(sum, 2);
%             carry = floor(sum / 2);
%         end
%         % Convert back to string
%         binaryStr = num2str(binaryVector);
%         binaryStr = strrep(binaryStr, ' ', ''); % Remove spaces
%     end
% end
% 
% % Example usage
% fraction = -10.625; % Example negative fractional number
% precision = 5;      % Number of bits for the fractional part
% binaryRepresentation = fractionalToBinary(fraction, precision);
% disp(['Binary representation: ', binaryRepresentation]);




% function binaryStr = fractionalToBinary(fraction, precision)
%     % Validate input
%     if fraction < 0
%         error('Input must be a non-negative fractional number.');
%     end
% 
%     % Initialize variables
%     integerPart = floor(fraction);
%     fractionalPart = fraction - integerPart;
%     binaryStr = '';
% 
%     % Convert integer part to binary
%     binaryStr = [binaryStr, dec2bin(integerPart), '.'];
% 
%     % Convert fractional part to binary
%     while precision > 0
%         fractionalPart = fractionalPart * 2;
%         bit = floor(fractionalPart);
%         binaryStr = [binaryStr, num2str(bit)];
%         fractionalPart = fractionalPart - bit;
%         precision = precision - 1;
%     end
% end
% 
% % Example usage
% fraction = 10.625; % Example fractional number
% precision = 5;     % Number of bits for the fractional part
% binaryRepresentation = fractionalToBinary(fraction, precision);
% disp(['Binary representation: ', binaryRepresentation]);




function binaryStr = fractionalToFixedPointBinary(fraction)
    % Define fixed-point parameters
    totalBits = 16;  % Total bits
    signBits = 1;    % Sign bit
    integerBits = 9; % Integer bits
    fractionalBits = 6; % Fractional bits

    % Check for NaN
    if isnan(fraction)
        binaryStr = 'NaN'; % Return 'NaN' as a string for NaN input
        return;
    end

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
        binaryVector = [str2num(binaryStr(:))'; 0]; % Append a zero for sign bit
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

% Example usage
fraction = -10.625; % Example negative fractional number
binaryRepresentation = fractionalToFixedPointBinary(fraction);
disp(['Binary representation: ', binaryRepresentation]);






