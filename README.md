# Fast_Hilbert_Transform
A pure Java library for fast computation of a 1D signal's hilbert transform.

## How it works
This library uses an implementation of the Fast Fourier Transform provided in the JTransfoms library in order to compute Hilbert transforms. This allows one to compute Hilbert transforms with signals of both even and odd lengths. The general algorithm is described [here](https://www.mathworks.com/help/signal/ref/hilbert.html).

