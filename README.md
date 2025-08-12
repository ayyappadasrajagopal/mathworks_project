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
    
