package com.quantimb;

import org.apache.commons.math3.transform.DftNormalization;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

public final class Hilbert {

  private final static DftNormalization dftNormalization = DftNormalization.STANDARD;
  private Hilbert() {}

  /**
   * Computes the hilbert transform of an input analytic signal using a Fast Fourier Transform.
   *
   * @param signal
   * @return
   */
  public static double[] computeHilbertTransform(double[] signal) {
    if (signal.length == 0) {
      return new double[0];
    }
    final DoubleFFT_1D transformer = new DoubleFFT_1D(signal.length);
    final double[] enlargedSignal = new double[signal.length * 2];

    for (int i = 0; i < signal.length; i++) {
      enlargedSignal[i] = signal[i];
    }

    transformer.realForwardFull(enlargedSignal);

    final double[] hVector = new double[enlargedSignal.length];
    hVector[0] = 1;
    final int middleElementIndex;
    if (enlargedSignal.length % 2 == 0) {
      hVector[enlargedSignal.length / 2] = 1;
      middleElementIndex = enlargedSignal.length / 2;
    } else {
      middleElementIndex = (enlargedSignal.length + 1) / 2;
    }
    for (int i = 1; i < middleElementIndex; i++) {
      hVector[i] = 2;
    }

    for (int i = 0; i < enlargedSignal.length; i++) {
      enlargedSignal[i] *= hVector[i];
    }

    transformer.complexInverse(enlargedSignal, true);

    return enlargedSignal;
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
