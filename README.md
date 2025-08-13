# mathworks_project


# 1. Fractional number to Binary 

### Step 1 – Understand the number

We have:
−10.625=−(10+0.625)
−10.625=−(10+0.625)
Step 2 – Choose a fixed-point format

Two’s complement for fractional numbers works just like for integers, but we fix how many integer bits and fractional bits we have.

Let’s say we pick an 8-bit fixed-point format:

    1 sign bit (in two’s complement, it’s part of the whole word)

    n integer bits (excluding sign bit)

    m fractional bits

For example: Q4.3 format →

    1 sign bit (part of 4 integer bits)

    3 fractional bits

So layout:
Sign/Integer bits	Fractional bits
4 bits	3 bits

Range for Q4.3:

    Min = –24−124−1 = –8 (Oops! This can’t hold –10)
    So we need more integer bits.

Let’s pick Q5.3:

    5 integer bits (including sign) → range –16 to +15.875

    3 fractional bits → resolution 1/8 = 0.125

That works for –10.625.
Step 3 – Write the positive version in binary

First, represent +10.625 in binary:
Integer part (10):
1010=10102
1010​=10102​

We need 5 integer bits (including sign bit for magnitude), so:
01010(for +10)
01010(for +10)
Fractional part (0.625):

Multiply repeatedly by 2:

    0.625×2=1.250.625×2=1.25 → 1

    0.25×2=0.50.25×2=0.5 → 0

    0.5×2=1.00.5×2=1.0 → 1

So fractional bits (3 bits) = 101.

Thus +10.625 in Q5.3 =
01010.1012
01010.1012​
Step 4 – Combine into fixed-point bits

Without the decimal point (just position it implicitly after 3 fractional bits):
01010 101(total 8 bits = Q5.3)
01010 101(total 8 bits = Q5.3)

That’s 01010101 in binary.
Step 5 – Convert to two’s complement for negative number

Two’s complement steps:

    Invert all bits:
    01010101→10101010
    01010101→10101010

    Add 1:
    10101010+1=10101011
    10101010+1=10101011

So –10.625 in Q5.3 two’s complement =
101010112
101010112​
Step 6 – Check

Interpret 10101011 as two’s complement Q5.3:

    Sign bit = 1 → negative

    Invert: 10101011 → 01010100

    Add 1: 01010100 + 1 = 01010101

    Place binary point after 3 bits: 01010.101 = 10.625

    Negate → –10.625

Final Answer:

    Q format: Q5.3 (8-bit total, 5 integer incl. sign, 3 fractional)

    Binary: 10101011

    Hex: AB



### Example with 

–10.625 for a 16-bit fixed-point format with:

    1 sign bit (part of the integer field in two’s complement)

    9 integer bits total (including the sign bit) → can represent integers from –512 to +511

    6 fractional bits → resolution = 1/64=0.0156251/64=0.015625

That means the format is Q9.6 (total 16 bits).
Step 1 – Positive version of the number

We first represent +10.625 in Q9.6.
Integer part:
1010=10102
1010​=10102​

In 9 integer bits (including sign for positive):
000001010(binary for +10)
000001010(binary for +10)
Fractional part (0.625):

Multiply by 2 repeatedly for 6 fractional bits:

    0.625×2=1.250.625×2=1.25 → 1

    0.25×2=0.50.25×2=0.5 → 0

    0.5×2=1.00.5×2=1.0 → 1

    0.0×2=0.00.0×2=0.0 → 0

    0.0×2=0.00.0×2=0.0 → 0

    0.0×2=0.00.0×2=0.0 → 0

Fractional bits = 101000.
Step 2 – Combine integer and fractional parts

+10.625 in Q9.6:
000001010.1010002
000001010.1010002​

Remove the binary point → 000001010101000 (16 bits total):
0000010101010000
0000010101010000

(We add a zero at the end because we have 9 integer bits + 6 fractional bits = 15 bits, plus the sign bit in position 15 already included.)

Let’s verify the placement carefully:

    Bits 15..7 → 9 integer bits (sign included)

    Bits 6..0 → 6 fractional bits

So:

Integer bits:   000001010   (binary 10)
Fractional bits:101000

Full 16-bit = 0000010101010000.
Step 3 – Two’s complement to make it negative

    Invert bits:

0000010101010000
→ 1111101010101111

Add 1:

    1111101010101111
    +                1
    = 1111101010110000

Step 4 – Final representation

    Binary: 1111101010110000

    Hex: 0xFA B0

    Q format: Q9.6

Step 5 – Check

Interpret 1111101010110000 as Q9.6:

    Sign bit = 1 → negative

    Invert: 0000010101001111

    Add 1: 0000010101010000

    Split: integer = 000001010 = 10, fraction = 101000 = 0.625

    → 10.625, negate → –10.625 ✅


    1) Top-level script flow

    You define two short messages:

    msg = ['Hello world 1'; 'Hello world 2'];

    MATLAB pads rows to equal length, so each row is a fixed-length string.

    You set radio knobs once:

        Platform: B200

        Serial/address: '30FE805'

        Gain: 35 dB

        Center frequency: 915 MHz

        Stop time (total): 100 s

        Sample rate (baseband Fs): 1 MHz

    You loop over each message, split the total time evenly between them, build a parameter struct for that message, create the modulator and radio, and transmit until the per-message stop time elapses. You also count hardware underruns (i.e., times the host didn’t feed samples fast enough).

2) sdruqpsktransmitter_init: builds all link parameters
Modulation & timing

    ModulationOrder = 4 → QPSK (2 bits/symbol).

    Interpolation = 2 → 2 samples/symbol at baseband after pulse shaping.

    Fs = sampleRate = 1e6 samples/s.

    Rsym = Fs / Interpolation = 1e6 / 2 = 500 ksym/s.

    Tsym = 1 / Rsym = 2 µs per symbol.

Implication: Raw (uncoded) bit rate = 2 bits/sym × 500 ksym/s = 1 Mb/s (before headers and framing).
Frame format (header + payload)

    Header: two Barker sequences back-to-back.

        Barker(13): [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]

        Header length = 13 × 2 = 26 symbols.

    Payload generation:

        You pass one base string (e.g., 'Hello world 1') for this run.

        Code builds 100 messages by appending a space + 3-digit counter + newline:

    "Hello world 1 000\n", "Hello world 1 001\n", … "Hello world 1 099\n"

    MessageLength = length(base) + 5 (space + 3 digits + newline).

    Characters are encoded in 7-bit ASCII: de2bi(..., 7, 'left-msb').

    PayloadLength = NumberOfMessage × MessageLength × 7 (bits).

FrameSize (in QPSK symbols) =

    (HeaderLength + PayloadLength) / log2(M)
    = (26 + PayloadLength) / 2

With your specific numbers (for intuition):

    If length('Hello world 1') = 13, then MessageLength = 18.

    Payload bits = 100 × 18 × 7 = 12,600 bits.

    FrameSize = (26 + 12,600)/2 = 6,313 symbols.

    Samples per frame (baseband) = Interpolation × FrameSize = 2 × 6,313 = 12,626 samples.

    Frame time = FrameSamples / Fs = 12,626 / 1e6 ≈ 12.626 ms.

Transmit signal shaping

    Root-Raised Cosine (RRC) shaping:

        RolloffFactor = 0.5

        RaisedCosineFilterSpan = 10 symbols (filter length ≈ 10 × Interpolation + 1 taps inside the object).

    Scrambler: base 2 with polynomial [1 1 1 0 1] (i.e., 1+D+D2+D41+D+D2+D4) and zero initial state—randomizes bit patterns to avoid long runs.

USRP clocking & rates

    MasterClockRate depends on platform:

        B200/B210 → 20 MHz

    The radio DAC/FPGA interface runs at:

        USRPFrontEndSampleRate = Rsym × 2 = 1e6 (equals your Fs).

        USRPInterpolationFactor = MasterClockRate / USRPFrontEndSampleRate = 20e6 / 1e6 = 20 (must be an integer—good).

    USRPFrameLength = Interpolation × FrameSize (samples per frame; matches the modulator’s output length).

    USRPFrameTime = USRPFrameLength / USRPFrontEndSampleRate (≈ 12.626 ms with the numbers above).

    StopTime is the per-message transmit duration (your total 100 s is split between the two messages by round(100/2) → 50 s each).

3) runSDRuQPSKTransmitter: creates objects and streams samples

    Uses persistent objects so they’re constructed once per run:

        QPSKTransmitter (MathWorks System object from the example set):

            Takes your bit stream, prepends Barker×2 header, scrambles, Gray-maps to QPSK, RRC-filters, and upsamples by Interpolation = 2, producing complex baseband at Fs = 1 MHz.

            Each hTx() call emits one frame (12,626 complex samples with your params).

        comm.SDRuTransmitter:

            Configured for B200 with SerialNum, MasterClockRate = 20 MHz, CenterFrequency = 915 MHz, Gain = 35 dB, InterpolationFactor = 20.

            Sends complex baseband to the radio; the USRP upsamples by the interpolation factor to the master clock.

    Main loop:

        data = hTx(); → generate the next frame.

        tunderrun = radio(data); → push frame to USRP; returns underrun count (0/1 or count).

        Accumulate underrun and advance currentTime by USRPFrameTime.

        Repeat until currentTime >= StopTime.

Underruns: If they’re non-zero, the host isn’t keeping up (USB/CPU scheduling). You’ll hear/see gaps over the air.
4) What to tweak (and what each knob does)

    Carrier frequency: USRPCenterFrequency (Hz). You used 915e6 (ISM band in many regions).

    Output level: USRPGain (dB). Increase carefully to avoid distortion/overdrive.

    User data rate:

        Baseband sample rate Fs controls Rsym = Fs/Interpolation.

        QPSK raw bit rate = 2 × Rsym.

        E.g., Fs=1e6, Interpolation=2 → Rsym=500 ksps → 1 Mb/s raw (before headers).

    Spectral shape: RolloffFactor and RaisedCosineFilterSpan.

    Frame content:

        Change SimParams.NumberOfMessage, SimParams.Message, or the format in sprintf to alter payload/length.

        Keep PayloadLength a multiple of 2 bits (QPSK) → the code already ensures that via the division by log2(M).

5) Common pitfalls & quick fixes

    No RF output / wrong device: Make sure SerialNum matches uhd_find_devices (or the USRP Radio Support Package “Find radios” tool). For X/N series, use IP instead.

    Non-integer interpolation: MasterClockRate / USRPFrontEndSampleRate must be an integer (it is: 20 MHz / 1 MHz = 20).

    Underruns:

        Try a slightly lower Fs or a smaller frame (reduce payload) to ease host load.

        Close other heavy apps/plots; disable debug prints inside the loop.

    Regulatory: 915 MHz usage depends on your region—ensure you’re compliant with local regulations, bandwidth, and power limits.

6) Quick mental “data budget” with your settings

    One frame ≈ 12.626 ms.

    Raw payload per frame ≈ FrameSize × 2 = ~12,626 bits (but includes header & 7-bit chars).

    At ~79 frames/s, you’re pushing ≈ 1 Mbit/s QPSK raw (pre-overheads).
    
