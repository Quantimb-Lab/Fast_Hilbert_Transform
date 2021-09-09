package com.quantimb;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

public final class Hilbert {

  private Hilbert() {}

  /**
   * <p>
   * Computes the hilbert transform of an input analytic signal using a Fast Fourier Transform.
   * </p>
   *
   * <p>
   * The outputted data has the form:
   * </p>
   *
   * <pre>
   * a[2*k] = Re[k],
   * a[2*k+1] = Im[k], 0<=k<n
   * </pre>
   *
   * @param signal
   * @return the hilbert transform, with data structured as specified above.
   */
  public static double[] computeHilbertTransform(double[] signal) {
    if (signal == null) {
      throw new NullPointerException("input signal must not be null");
    }
    final double[] hilbertTransform;
    if (signal.length == 0) {
      hilbertTransform = new double[0];
    }
    else {
      final DoubleFFT_1D transformer = new DoubleFFT_1D(signal.length);
      final double[] fourierSignal = new double[signal.length * 2];

      for (int i = 0; i < signal.length; i++) {
        fourierSignal[i] = signal[i];
      }

      transformer.realForwardFull(fourierSignal);

      final double[] hVector = createHVector(fourierSignal);

      for (int i = 0; i < fourierSignal.length; i++) {
        fourierSignal[i] *= hVector[i];
      }

      transformer.complexInverse(fourierSignal, true);
      hilbertTransform = fourierSignal;
    }

    return hilbertTransform;
  }

  /**
   * <p>
   * Creates the H vector for Hilbert transform
   * </p>
   *
   * <p>
   * If we seek to send across a list of complex numbers [(a+bi), (c+di), ...], then the input
   * fourier signal should have the following form: [a,b,c,d,...]
   * </p>
   *
   * @param fourierSignal
   * @return
   */
  private static double[] createHVector(double[] fourierSignal) {
    final double[] hVector = new double[fourierSignal.length];
    setPairOfDoubleValues(hVector, 0, 1);

    final int middleComplexElementIndex;
    final int numComplexElements = fourierSignal.length / 2;

    if (numComplexElements % 2 == 0) {
      middleComplexElementIndex = (numComplexElements) / 2;
      setPairOfDoubleValues(hVector, middleComplexElementIndex, 1);
    } else {
      middleComplexElementIndex = (numComplexElements + 1) / 2;
    }

    for (int i = 1; i < middleComplexElementIndex; i++) {
      setPairOfDoubleValues(hVector, i, 2);
    }
    return hVector;
  }

  /**
   * Assumes arr is an array of complex numbers [(a+bi), (c+di), ...] structured as [a,b,c,d,...]
   *
   * @param arr
   * @param index
   * @param val
   */
  private static void setPairOfDoubleValues(double[] arr, int index, int val) {
    arr[2 * index] = val;
    arr[2 * index + 1] = val;
  }

  public static double[] computeSignalEnvelope(double[] hilbertTransform) {
    final double[] envelope = new double[hilbertTransform.length / 2];
    for (int i = 0; i < envelope.length; i++) {
      envelope[i] = Math
          .sqrt(Math.pow(hilbertTransform[2 * i], 2) + Math.pow(hilbertTransform[2 * i + 1], 2));
    }
    return envelope;
  }
}
