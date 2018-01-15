# sinc-errorbound-infinite
Sinc approximation with its rigorous error bound over infinite interval

## Overview
These programs approximate the following four functions:

(1) f(t) = sqrt(1 + tanh^2(arcsinh t)) / (1 + t^2)  
(2) f(t) = sqrt(t) sqrt(1 + tanh^2(log t)) / (1 + t^2)  
(3) f(t) = t^(pi/4) exp(-t)  
(4) f(t) = sqrt(cos(3 arcsinh t) + cosh(pi)) / (1 + t^2)

Approximation is done by means of the Sinc approximation combined with
a proper variable transformation. Programs that begin with SE_Sinc use
variable transformations given by Stenger [1], whereas other programs
use the so-called double-exponential transformation described in [2, 3].

Each program investigates maximum approximation error among selected
sampling points increasing N as N = 2, 7, 12, 17, 22, ..., then outputs
N, maximum error, and its error bound. Note that error bounds by
OLD_DE_Sinc_3.c and DE_Sinc_4.c are not correct (just for reference).

## Results
Outputs by those programs are stored in data/ directory, with .dat extension.
Gnuplot programs and created eps graphs are also stored in the directory.

computation environment:

OS: Mac OS X 10.9  
CPU: 1.7 GHz Intel Core i7  
Memory: 8 GB 1600 MHz DDR3  
Compiler: Apple LLVM version 6.0  

## References
1. F. Stenger:
 Numerical Methods Based on Sinc and Analytic Functions, Springer-Verlag,
 New York, 1993.
2. K. Tanaka, M. Sugihara, K. Murota:
 Function classes for successful DE-Sinc approximations, Mathematics of
 Computation, Vol. 78 (2009), pp. 1553--1571.
3. T. Okayama:
 Error estimates with explicit constants for the Sinc approximation over
 infinite intervals, Applied Mathematics and Computation, Vol. 319 (2018),
 pp. 125--137.
