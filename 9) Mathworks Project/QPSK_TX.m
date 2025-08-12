%% % mESSAGE TRANSMISSION USING QPSK with USRP Hardware

clc;
close all;
clear all;

% contnous message transmission of predefined set of messages
msg = ['Hello world 1';'Hello world 2'];

platform = "B200"; % NI USRP device
address  = '30FE805';  % NI USRP device address
USRPGain            = 35;  % Set USRP radio gain
USRPCenterFrequency = 915000000;  % Set USRP radio center frequency
stopTimeb           = 100;  % USRP radio transmit time in seconds
sampleRate          = 1000000;  % Sample rate of transmitted signal

for i = 1:size(msg,1)
stopTime = round(stopTimeb/size(msg,1)); 
message = msg(i,:);
% Transmitter parameter structure
prmQPSKTransmitter = sdruqpsktransmitter_init(platform, address, sampleRate, USRPCenterFrequency, ...
    USRPGain, stopTime,message);
underruns = runSDRuQPSKTransmitter(prmQPSKTransmitter);
fprintf('Total number of underruns = %d.\n', underruns);
end

function SimParams = sdruqpsktransmitter_init(platform,address,sampleRate,....
    centerFreq,gain,captureTime,message)
%   Copyright 2012-2023 The MathWorks, Inc.

%% General simulation parameters
SimParams.ModulationOrder = 4; % QPSK alphabet size
SimParams.Interpolation   = 2; % Interpolation factor
SimParams.Decimation      = 1; % Decimation factor
SimParams.Fs              = sampleRate; % Sample rate
SimParams.Rsym            = sampleRate/SimParams.Interpolation; % Symbol rate in Hertz
SimParams.Tsym            = 1/SimParams.Rsym; % Symbol time in sec

%% Frame Specifications
% [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ...];
SimParams.BarkerCode      = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];     % Bipolar Barker Code
SimParams.BarkerLength    = length(SimParams.BarkerCode);
SimParams.HeaderLength    = SimParams.BarkerLength * 2;                   % Duplicate 2 Barker codes to be as a header
SimParams.Message         = message;
SimParams.MessageLength   = length(SimParams.Message) + 5;                % 'Hello world 000\n'...
SimParams.NumberOfMessage = 100;                                          % Number of messages in a frame
SimParams.PayloadLength   = SimParams.NumberOfMessage * SimParams.MessageLength * 7; % 7 bits per characters
SimParams.FrameSize       = (SimParams.HeaderLength + SimParams.PayloadLength) ...
    / log2(SimParams.ModulationOrder);                                    % Frame size in symbols
SimParams.FrameTime       = SimParams.Tsym*SimParams.FrameSize;

%% Tx parameters
SimParams.RolloffFactor     = 0.5;                                        % Rolloff Factor of Raised Cosine Filter
SimParams.ScramblerBase     = 2;
SimParams.ScramblerPolynomial           = [1 1 1 0 1];
SimParams.ScramblerInitialConditions    = [0 0 0 0];
SimParams.RaisedCosineFilterSpan = 10; % Filter span of Raised Cosine Tx Rx filters (in symbols)

%% Message generation
msgSet = zeros(100 * SimParams.MessageLength, 1); 
for msgCnt = 0 : 99
    msgSet(msgCnt * SimParams.MessageLength + (1 : SimParams.MessageLength)) = ...
        sprintf('%s %03d\n', SimParams.Message, msgCnt);
end
bits = de2bi(msgSet, 7, 'left-msb')';
SimParams.MessageBits = bits(:);

%% USRP transmitter parameters
SimParams.Platform                      = platform;
SimParams.Address                       = address;

switch platform
  case {'B200','B210'}
    SimParams.MasterClockRate = 20e6;           % Hz
  case {'X300','X310'}
    SimParams.MasterClockRate = 200e6;          % Hz
  case {'N300','N310'}
    SimParams.MasterClockRate = 125e6;          % Hz
  case {'N320/N321'}
    SimParams.MasterClockRate = 200e6;          % Hz
  case {'N200/N210/USRP2'}
    SimParams.MasterClockRate = 100e6;          % Hz
  otherwise
    error(message('sdru:examples:UnsupportedPlatform', ...
      platform))
end
SimParams.USRPCenterFrequency       = centerFreq;
SimParams.USRPGain                  = gain;
SimParams.USRPFrontEndSampleRate    = SimParams.Rsym * 2; % Nyquist sampling theorem
SimParams.USRPInterpolationFactor   = SimParams.MasterClockRate/SimParams.USRPFrontEndSampleRate;
SimParams.USRPFrameLength           = SimParams.Interpolation * SimParams.FrameSize;

% Experiment Parameters
SimParams.USRPFrameTime = SimParams.USRPFrameLength/SimParams.USRPFrontEndSampleRate;
SimParams.StopTime = captureTime;
end

function underrun = runSDRuQPSKTransmitter(prmQPSKTransmitter)
%#codegen

%   Copyright 2012-2023 The MathWorks, Inc.

    persistent hTx radio
    if isempty(hTx)
        % Initialize the components
        % Create and configure the transmitter System object
        hTx = QPSKTransmitter(...
            'UpsamplingFactor',             prmQPSKTransmitter.Interpolation, ...
            'RolloffFactor',                prmQPSKTransmitter.RolloffFactor, ...
            'RaisedCosineFilterSpan',       prmQPSKTransmitter.RaisedCosineFilterSpan, ...
            'MessageBits',                  prmQPSKTransmitter.MessageBits, ...
            'MessageLength',                prmQPSKTransmitter.MessageLength, ...
            'NumberOfMessage',              prmQPSKTransmitter.NumberOfMessage, ...
            'ScramblerBase',                prmQPSKTransmitter.ScramblerBase, ...
            'ScramblerPolynomial',          prmQPSKTransmitter.ScramblerPolynomial, ...
            'ScramblerInitialConditions',   prmQPSKTransmitter.ScramblerInitialConditions);

        % Create and configure the SDRu System object. Set the SerialNum for B2xx
        % radios and IPAddress for X3xx, N2xx, and USRP2 radios. MasterClockRate
        % is not configurable for N2xx and USRP2 radios.
        switch prmQPSKTransmitter.Platform
            case {'B200','B210'}
                radio = comm.SDRuTransmitter(...
                    'Platform',             prmQPSKTransmitter.Platform, ...
                    'SerialNum',            prmQPSKTransmitter.Address, ...
                    'MasterClockRate',      prmQPSKTransmitter.MasterClockRate, ...
                    'CenterFrequency',      prmQPSKTransmitter.USRPCenterFrequency, ...
                    'Gain',                 prmQPSKTransmitter.USRPGain, ...
                    'InterpolationFactor',  prmQPSKTransmitter.USRPInterpolationFactor);
            case {'X300','X310'}
                radio = comm.SDRuTransmitter(...
                    'Platform',             prmQPSKTransmitter.Platform, ...
                    'IPAddress',            prmQPSKTransmitter.Address, ...
                    'MasterClockRate',      prmQPSKTransmitter.MasterClockRate, ...
                    'CenterFrequency',      prmQPSKTransmitter.USRPCenterFrequency, ...
                    'Gain',                 prmQPSKTransmitter.USRPGain, ...
                    'InterpolationFactor',  prmQPSKTransmitter.USRPInterpolationFactor);
            case {'N200/N210/USRP2'}
                radio = comm.SDRuTransmitter(...
                    'Platform',             prmQPSKTransmitter.Platform, ...
                    'IPAddress',            prmQPSKTransmitter.Address, ...
                    'CenterFrequency',      prmQPSKTransmitter.USRPCenterFrequency, ...
                    'Gain',                 prmQPSKTransmitter.USRPGain, ...
                    'InterpolationFactor',  prmQPSKTransmitter.USRPInterpolationFactor);
            case {'N300','N310'}
                radio = comm.SDRuTransmitter(...
                    'Platform',             prmQPSKTransmitter.Platform, ...
                    'IPAddress',            prmQPSKTransmitter.Address, ...
                    'MasterClockRate',      prmQPSKTransmitter.MasterClockRate, ...
                    'CenterFrequency',      prmQPSKTransmitter.USRPCenterFrequency, ...
                    'Gain',                 prmQPSKTransmitter.USRPGain, ...
                    'InterpolationFactor',  prmQPSKTransmitter.USRPInterpolationFactor);
            case {'N320/N321'}
                radio = comm.SDRuTransmitter(...
                    'Platform',             prmQPSKTransmitter.Platform, ...
                    'IPAddress',            prmQPSKTransmitter.Address, ...
                    'MasterClockRate',      prmQPSKTransmitter.MasterClockRate, ...
                    'CenterFrequency',      prmQPSKTransmitter.USRPCenterFrequency, ...
                    'Gain',                 prmQPSKTransmitter.USRPGain, ...
                    'InterpolationFactor',  prmQPSKTransmitter.USRPInterpolationFactor);                
        end
    end    
    
    cleanupTx = onCleanup(@()release(hTx));
    cleanupRadio = onCleanup(@()release(radio));

    currentTime = 0;
    underrun = uint32(0);

    %Transmission Process
    while currentTime < prmQPSKTransmitter.StopTime
        % Bit generation, modulation and transmission filtering
        data = hTx();

        % Data transmission
        tunderrun = radio(data);
        underrun = underrun + tunderrun;

        % Update simulation time
        currentTime=currentTime+prmQPSKTransmitter.USRPFrameTime;
    end

end
