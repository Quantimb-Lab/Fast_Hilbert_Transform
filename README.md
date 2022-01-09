# Fast_Hilbert_Transform
A pure Java library for fast computation of a 1D signal's hilbert transform.

## How it works
This library uses an implementation of the Fast Fourier Transform provided in the JTransfoms library in order to compute Hilbert transforms. This allows one to compute Hilbert transforms with signals of both even and odd lengths. The general algorithm is described [here](https://www.mathworks.com/help/signal/ref/hilbert.html). This allows for computation of a Hilbert transform in O(n log n), where n is a signal's length

## How to use it
```Java
# represent a 1D signal as an array of doubles
final double[] signal = {.....};

# compute the Hilbert transform of the signal
final double[] hilbertTransform = Hilbert.computeHilbertTransform(signal);

# compute the amplitude envelope of the Hilbert transformed signal
final double[] amplitudeEnvelope = Hilbert.computeSignalEnvelope(hilbertTransform);
```

## Importing into a maven project
In your pom file, add the following plugin:

```
<plugin>
    <groupId>io.quantimb_lab.quantimb</groupId>
    <artifactId>hilbert_transform</artifactId>
    <version>1.0-SNAPSHOT</version>
</plugin>
```