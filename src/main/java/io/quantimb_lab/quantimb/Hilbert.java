  /**
  * Copyright 2021 Quantimb Lab
  *
  * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
  * 
  * You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
  *
  * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  *
  * See the License for the specific language governing permissions and limitations under the License.
  *
  */
  
package io.quantimb_lab.quantimb;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

  /**
   * <h1>Hilbert Tansform (FFT Implementation)</h1>
   * 
   * <p>
   * The Hilbert Class is developed to compute the hilbert transform of an input signal using a Fast Fourier Transform.
   * The Amplitude Envelope can then be found using the analytical signal created by the Hilbert Transform of the input signal.
   * </p>
   *
   * @author  Yaniv Khaslavsky, Reza Yaghoubi
   * @version 1.0
   */
   
public final class Hilbert {

  private Hilbert() {}

  /**
   * 
   * <p>
   * Computes the hilbert transform of an input signal using a Fast Fourier Transform.
   * It produces an analytical signal from which the amplitude envelope can then be found.
   * </p>
   *
   * <p>
   * The outputted data has the form:
   * </p>
   *
   * <pre>
   * a[2*k] = Re[k],
   * a[2*k+1] = Im[k], 0&lt;=k&lt;n
   * </pre>
   * @param signal
   * @return the hilbert transform, with data structured as specified above.
   *
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
   * @return hVector
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


  /**
   * This method computes the Amplitude Envelope of the input signal for which the Hilbert Transform is already computed
   *
   * @param hilbertTransform (the hilbert Transform of the input signal)/the analytical signal
   * @return Amplitude Envelope
   * 
   */
  public static double[] computeSignalEnvelope(double[] hilbertTransform) {
    final double[] envelope = new double[hilbertTransform.length / 2];
    for (int i = 0; i < envelope.length; i++) {
      envelope[i] = Math
          .sqrt(Math.pow(hilbertTransform[2 * i], 2) + Math.pow(hilbertTransform[2 * i + 1], 2));
    }
    return envelope;
  }
}
