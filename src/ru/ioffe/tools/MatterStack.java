/**
 * 
 */
package ru.ioffe.tools;

import java.util.ArrayList;

import javax.swing.JFrame;

import org.opensourcephysics.frames.PlotFrame;

import ru.ioffe.semiconductor.Utils;

/**
 * @author Evgeny Shevchenko
 * 
 */
public strictfp class MatterStack {

	private double ENERGY_STEP = Utils.ev2erg(0.001); // 0.0001 eV in ergs
	private int X_DOTS_DENSITY = 40; // per 1 nm
	private double[] eigenEnergy = null;
	private double[] energyRange = null;
	private ArrayList<MatterInterval> intervals = new ArrayList<MatterInterval>();
	private double y_shift = 0;
	// private boolean fictitiousRegionReady = false;

	// private boolean VERBOSE = false;

	private boolean VERBOSE = true;

	/**
	 * 
	 */
	public MatterStack() {
	}

	/**
	 * 
	 * @param _intervals
	 */
	public MatterStack(ArrayList<MatterInterval> _intervals) {

		this.intervals = _intervals;
		double width = 0;
		y_shift = Utils.ev2erg(intervals.get(0).getY());

		for (MatterInterval inl : intervals) {
			inl.setXStart(width);
			width += inl.getWidth();
			inl.setY(Utils.ev2erg(inl.getY()));
			if (y_shift > inl.getY()) {
				y_shift = inl.getY();
			}
		}

		for (MatterInterval inl : intervals) {
			inl.setY(inl.getY() - y_shift);
		}

		MatterInterval last = intervals.get(intervals.size() - 1);
		// intervals.add(new MatterInterval(width, last.getWidth() + width,
		// last.getY(), last.getMass(), last.getDiel()));
		intervals.add(new MatterInterval(width, width, last.getY(), last
				.getMass(), last.getDiel()));
		intervals.get(intervals.size() - 1).setLabel("Fictitious");
	}

	@Override
	public String toString() {
		StringBuffer s = new StringBuffer();
		int counter = 0;
		for (MatterInterval mi : intervals) {
			s.append("" + (counter++) + ": ");
			s.append(mi.toString() + "\n");
		}
		return s.toString();
	}

	public MatterInterval getInterval(int index) {
		return intervals.get(index);
	}

	public MatterInterval getFictitiousInterval() {
		return intervals.get(intervals.size() - 1);
	}

	/**
	 * Takes index of preceding interval (between index and index+1) and energy
	 * value
	 * 
	 * @return ComplexMatrix
	 * @throws Exception
	 */
	public ComplexMatrix getLocalPropagationMatrix(int index, double energy) {
		MatterInterval inl = intervals.get(index);
		MatterInterval next = intervals.get(index + 1);
		Complex p;
		// Double.POSITIVE_INFINITY is handled here and it assigns to real part
		// of Complex if appears
		p = Complex.sqrt((energy - next.getY()) / (energy - inl.getY())
				* inl.getMass() / next.getMass());

		Complex kcp = p.add(1); // 1 + p
		Complex kcm = p.reverseSign().add(1); // 1 - p
		Complex kd = next.getWaveNumber(energy).times(next.getWidth());

		if (VERBOSE) {
			System.out.println("p   = " + p);
			System.out.println("kcp = " + kcp);
			System.out.println("kcm = " + kcm);
			System.out.println("kd  = " + kd);
		}

		Complex ipos = new Complex(0, 1);
		Complex ineg = new Complex(0, -1);
		try {
			return new ComplexMatrix(new Complex[][] {
					{ kcp.times(Complex.exp(ipos.times(kd))),
							kcm.times(Complex.exp(ipos.times(kd))) },
					{ kcm.times(Complex.exp(ineg.times(kd))),
							kcp.times(Complex.exp(ineg.times(kd))) } })
					.multiplyEach(0.5);
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Get array of complex coefficients of exponential 1-dimensional wave
	 * function
	 * 
	 * @param index
	 * @param energy
	 * @return
	 */
	public ComplexMatrix getPropagationMatrix(int index, double energy) {
		ComplexMatrix mx = ComplexMatrix.identityMatrix(2);
		for (int i = 0; i < index; i++) {
			try {
				mx = getLocalPropagationMatrix(i, energy).multiply(mx);
				if (VERBOSE) {
					System.out.println("> > > In getPropagationMatrix(...)");
					System.out.println("Matrix on index " + index + ", energy "
							+ energy + ", interation " + i);
					System.out.println(mx);
					System.out.println("Current propagation matrix = "
							+ getLocalPropagationMatrix(i, energy));
				}
			} catch (NullPointerException e) {
				System.out
						.println("Null received in ComplexMatrix multiplication in getPropagationMatrix(...)");
				e.printStackTrace();
			} catch (Exception e) {
				System.out
						.println("Something wrong on ComplexMatrix multiplication in getPropagationMatrix(...)");
				e.printStackTrace();
			}
		}
		return mx;
	}

	/**
	 * Get double array of energy range for search of eigen values. The range
	 * step is taken from ENERGY_STEP
	 * 
	 * @return
	 */
	public double[] getEnergyRange() {
		// EnergyRange is stored in instance in order not to generate it
		// multiple times
		if (energyRange != null) {
			return energyRange;
		}

		double max = 0;
		double max2 = 0;

		for (MatterInterval mi : intervals) {
			if (max <= mi.getY()) {
				max2 = max;
				max = mi.getY();
			}
		}
		int rangeSize = (int) (max2 / ENERGY_STEP);
		if (VERBOSE) {
			System.out.println("Energy max:  " + Utils.erg2ev(max));
			System.out.println("Energy max2: " + Utils.erg2ev(max2));
			System.out.println("Energy step: " + ENERGY_STEP);
			System.out.println("Energy range size: " + rangeSize);
		}
		energyRange = new double[rangeSize];

		for (int i = 0; i < rangeSize; i++) {
			energyRange[i] = i * ENERGY_STEP;
		}
		return energyRange;
	}

	/**
	 * Value of propagation curve. Propagation maximums correspond to eigen
	 * states
	 * 
	 * @param index
	 *            index of a region/interval
	 * @param energy
	 *            energy in ergs
	 * @return
	 */
	// public double getPropagationValue(int index, double energy) {
	// if (VERBOSE) {
	// System.out.println("Propagation value for index" + index + " and energy "
	// + energy);
	// // System.out.println("" + getPropagationMatrix(index,
	// // energy)[0].abs());
	// // System.out.println("" + getPropagationMatrix(0,
	// // energy)[0].abs());
	// }
	// return Math.pow(getWavefunctionCoeffs(index, energy)[0].abs() /
	// getWavefunctionCoeffs(0, energy)[0].abs(), 2);
	// }
	public double getPropagationValue(double energy) {
		if (VERBOSE) {
			System.out.println("Propagation value for energy " + energy);
		}
		return getPropagationMatrix(intervals.size() - 1, energy).getElem(1, 1)
				.abs();
	}

	/**
	 * Derivative of the propagation curve
	 * 
	 * @param index
	 *            index of a region/interval
	 * @param energy
	 *            energy in ergs
	 * @return
	 */
	public double getPropagationValueDerivative(double energy) {
		return (getPropagationValue(energy + Function.CALC_TOLERANCE) - getPropagationValue(energy
				- Function.CALC_TOLERANCE))
				/ 2 / Function.CALC_TOLERANCE;
	}

	/*
	 * Needed only for Newton-Raphson method
	 * 
	 * public double getPropagationValueDerivative2(int index, double energy) {
	 * return (getPropagationValueDerivative(index, energy +
	 * Function.CALC_TOLERANCE) - getPropagationValueDerivative( index, energy -
	 * Function.CALC_TOLERANCE)) / 2 / Function.CALC_TOLERANCE; }
	 */

	/**
	 * Log of the value of propagation curve. Propagation maximums correspond to
	 * eigen states
	 * 
	 * @param index
	 *            index of a region/interval
	 * @param energy
	 *            energy in ergs
	 * @return
	 */
	private double getPropagationValueLog(double energy) {
		return Math.log(getPropagationValue(energy));
	}

	/**
	 * Derivative of the Log of propagation curve
	 * 
	 * @param index
	 *            index of a region/interval
	 * @param energy
	 *            energy in ergs
	 * @return
	 */
	@SuppressWarnings("unused")
	private double getPropagationValueLogDerivative(double energy) {
		return (getPropagationValueLog(energy + Function.CALC_TOLERANCE) - getPropagationValueLog(energy
				- Function.CALC_TOLERANCE))
				/ 2 / Function.CALC_TOLERANCE;
	}

	/**
	 * Get propagation energy in an interval between two energy values by
	 * bisection method
	 * 
	 * @param index
	 *            index of region/interval
	 * @param energy_min
	 *            minimal energy in ergs for the range to search in
	 * @param energy_max
	 *            maximum energy in ergs for the range to search in
	 * @return energy value in ergs if maximum is found
	 */
	public double getPropagationEnergyRoot(double energy_min, double energy_max) {
		System.out.println("ENTERING BISECTION in getPropagationEnergyRoot...");
		double energy = (energy_min + energy_max) / 2;
		// double dprop = getPropagationValueLogDerivative(index, energy);
		double dprop = getPropagationValueDerivative(energy);
		while (Math.abs(dprop) > Function.CALC_TOLERANCE) {
			System.out.println("BSF: " + dprop + "\tenergy_min: " + energy_min
					+ "\tenergy: " + energy + "\tenergy_max: " + energy_max);
			if (Math.abs(energy / (energy_min + energy_max) - 0.5) < Function.CALC_TOLERANCE) {
				return energy;
			}
			// if (dprop * getPropagationValueLogDerivative(index, energy_min) >
			// 0) {
			if (dprop * getPropagationValueDerivative(energy_min) > 0) {
				energy_min = energy;
				// } else if (dprop * getPropagationValueLogDerivative(index,
				// energy_max) > 0) {
			} else if (dprop * getPropagationValueDerivative(energy_max) > 0) {
				energy_max = energy;
			}
			if (energy == (energy_min + energy_max) / 2) {
				return energy;
			}
			System.out.println("Energy difference = "
					+ (energy - (energy_min + energy_max) / 2));
			energy = (energy_min + energy_max) / 2;
			// dprop = getPropagationValueLogDerivative(index, energy);
			dprop = getPropagationValueDerivative(energy);
		}
		return energy;
	}

	/**
	 * Get eigen energies for an interval
	 * 
	 * @param index
	 *            of an interval
	 * @return double[] of eigen energies in ergs
	 */
	// public double[] findIntervalEigenEnergies(int index) {
	// ArrayList<Double> energyList = new ArrayList<Double>();
	// double[] eRange = getEnergyRange();
	// for (int i = 1; i < eRange.length - 1; i++) {
	// // if (VERBOSE) {
	// // System.out.println("Index of interval: " + index + "    energy: "
	// // + eRange[i] + " erg, "
	// // + Utils.erg2ev(eRange[i]) + " eV");
	// // System.out.println("PropagationValueDerivative(" + (index - 1) +
	// // ", " + eRange[i] + ") = "
	// // + getPropagationValueDerivative(index - 1, eRange[i]));
	// // System.out.println("PropagationValueDerivative(" + (index) + ", "
	// // + eRange[i] + ") = "
	// // + getPropagationValueDerivative(index, eRange[i]));
	// // System.out.println("PropagationValueDerivative(" + (index + 1) +
	// // ", " + eRange[i] + ") = "
	// // + getPropagationValueDerivative(index + 1, eRange[i]));
	// // }
	// if (Double.isNaN(getPropagationValueDerivative(index, eRange[i - 1]))
	// || Double.isNaN(getPropagationValueDerivative(index, eRange[i]))
	// || Double.isNaN(getPropagationValueDerivative(index, eRange[i + 1]))) {
	// System.out.println("NaNs detected!");
	// continue;
	// } else if (getPropagationValue(index, eRange[i]) >
	// getPropagationValue(index, eRange[i - 1])
	// && getPropagationValue(index, eRange[i]) > getPropagationValue(index,
	// eRange[i + 1])
	// && getPropagationValueDerivative(index, eRange[i - 1]) *
	// getPropagationValueDerivative(index, eRange[i + 1]) < 0) {
	// // energyList.add(eRange[i]);
	// energyList.add(getPropagationEnergy(index, eRange[i - 1], eRange[i +
	// 1]));
	// // System.out.println("FOUND " + eRange[i] + "\tCALCULATED "
	// // + getPropagationEnergy(index, eRange[i - 1], eRange[i + 1]));
	// }
	// }
	// double[] energyArray = new double[energyList.size()];
	// for (int i = 0; i < energyList.size(); i++) {
	// energyArray[i] = Double.valueOf(energyList.get(i));
	// }
	// return energyArray;
	// }

	/**
	 * Build a plot of quantum wells propagation curves
	 */
	public void plotPropagationCurves() {
		double[] eRange = getEnergyRange();
		// X label, Y label, plot title
		PlotFrame frame = new PlotFrame("Energy, eV", 
										"Propagation, arb. units", 
										"QW Propagation Curves");
		// Filling the curve data
		for (int i = 0; i < eRange.length; i++) {
			frame.append(0, Utils.erg2ev(eRange[i]), getPropagationValueLog(eRange[i]));
		}
		// No markers
		frame.setMarkerShape(0, 0);

		// Window size 600 x 600
		frame.setSize(600, 600);

		// Dots are connected with a line
		frame.setConnected(true);

		// Materializing the plot
		frame.setVisible(true);

		// Behavior of the close button
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	/**
	 * Get eigen energies of the whole system
	 * 
	 * @param index
	 * @return double[] of eigen energies in ergs
	 */
	public double[] getEigenEnergy() {
		if (eigenEnergy != null) {
			return eigenEnergy;
		}
		ArrayList<Double> eigenEnergyList = new ArrayList<Double>();
		for (int i = 1; i < intervals.size() - 1; i++) {
			if (intervals.get(i).getY() <= intervals.get(i - 1).getY()
					&& intervals.get(i).getY() < intervals.get(i + 1).getY()) {
				// for (double e : findIntervalEigenEnergies(i)) {
				// if (!eigenEnergyList.contains(Double.valueOf(e))) {
				// eigenEnergyList.add(new Double(e));
				// }
				// }
			}
		}
		eigenEnergy = new double[eigenEnergyList.size()];
		for (int i = 0; i < eigenEnergy.length; i++) {
			eigenEnergy[i] = Double.valueOf(eigenEnergyList.get(i));
		}
		return eigenEnergy;
	}

	/**
	 * Print a column of eigen energy values of the stack
	 * 
	 * @param number
	 *            number of states starting with the lowest one (not implemented
	 *            yet)
	 */
	public void printEigenEnergy(int number) {
		double[] ee = getEigenEnergy();
		if (ee.length == 0) {
			System.out.println("No eigen states found :-(");
			return;
		}
		int count = Math.min(number, ee.length);
		for (int i = 0; i < count; i++) {
			System.out.println(Utils.erg2ev(ee[i]) * 1e3 + " meV");
		}
	}

	/**
	 * Get the whole MatterStack width starting from index=1 (thus excluding the
	 * fictitious interval) to the last one
	 * 
	 * @return
	 */
	public double getWidth() {
		if (intervals == null || intervals.size() == 0) {
			return 0.0;
		}
		return getSubStackWidth(1, intervals.size() - 1);
	}

	/**
	 * Get width of specified interval. Note that indexes start with 1 due to
	 * fictitious interval (0th)
	 * 
	 * @param index
	 * @return the interval width
	 */
	public double getIntervalWidth(int index) {
		return intervals.get(index).getWidth();
	}

	/**
	 * Get width of a part of the stack including both intervals
	 * 
	 * @param index1
	 *            number of starting interval
	 * @param index2
	 *            number of ending interval (included)
	 * @return common width of the sub-stack
	 */
	public double getSubStackWidth(int index1, int index2) {
		double width = 0;
		for (int i = index1; i <= index2; i++) {
			width += intervals.get(i).getWidth();
		}
		return width;
	}

	/**
	 * Get wave function coefficients for chosen region and energy
	 * 
	 * @param index
	 * @param energy
	 * @return
	 */
	public Complex[] getWavefunctionCoeffs(int index, double energy) {
		ComplexMatrix mx = getPropagationMatrix(index, energy);
		if (index == 1) {
			return new Complex[] { Complex.fromDouble(0), mx.getElem(1, 0) };
		}
		return new Complex[] { mx.getElem(0, 0), mx.getElem(1, 0) };
	}

	/**
	 * Get wave function for a state
	 * 
	 * @param energy
	 * @return double[] array of a wave function values
	 */
	public double[][] getWaveFunction(double energy) {

		// Number of dots in the whole MatterStack
		// FIXME magic number
		int numberOfDots = (int) (getWidth() * X_DOTS_DENSITY * 1000000000L);

		// TODO try-catch or think on long-or-int problem
		double[] waveFunction = new double[numberOfDots];
		double[] xAxis = new double[numberOfDots];
		double[] xLocal;

		int counter = 0;

		Complex[] coef; // {A, B}
		Complex pi = new Complex(0, 1);
		Complex mi = new Complex(0, -1);

		// Going through all intervals to build the wave funtion
		for (int k = 0; k < intervals.size(); k++) {

			MatterInterval inl = intervals.get(k);

			// FIXME magic number
			xLocal = inl.getGrid(X_DOTS_DENSITY * 1000000000L);

			for (int i = 0; i < xLocal.length; i++) {

				coef = getWavefunctionCoeffs(k, energy);

				// WF = | A * exp(i*k*x) + B * exp(-i*k*x) | -- doubles not
				// complex!
				waveFunction[i + counter] = (coef[0].times(Complex.exp(pi
						.times(inl.getWaveNumber(energy).times(
								xLocal[i] - inl.getXStart()))))).add(
						coef[1].times(Complex.exp(mi.times(inl.getWaveNumber(
								energy).times(xLocal[i] - inl.getXStart())))))
						.abs();

				xAxis[i + counter] = xLocal[i];

			}
			counter += xLocal.length;
		}
		return new double[][] { xAxis, waveFunction };
	}

	/**
	 * Get potential profile
	 * 
	 * @return double[][] of the potential's x and y
	 */
	public double[][] getPotentialProfile() {
		double[] x = new double[intervals.size()];
		double[] y = new double[x.length];
		for (int i = 0; i < y.length; i++) {
			x[i] = intervals.get(i).getX();
			y[i] = intervals.get(i).getY();
		}
		return new double[][] { x, y };
	}

	// private double[] reduceWaveFunction() {}

	/**
	 * @param args
	 */
	@SuppressWarnings("unused")
	public static void main(String[] args) {

		double Eg = 4.0;
		double dEg = 0.1;
		double qww = 2e-7;
		double qwb = qww / 2;
		double b1w = qww * 1e2;
		double b2w = b1w;
		double diel = 10;
		MatterStack a;

		// a.appendNewMatterInterval(b1w, Eg, 0.20, 9.9);
		// a.appendNewMatterInterval(qww, Eg - dEg, 0.25, 10);
		// a.appendNewMatterInterval(qww, Eg - dEg / 2, 0.25, 10);
		// a.appendNewMatterInterval(qwb, Eg, 0.20, 10);
		// a.appendNewMatterInterval(b1w, Eg, 0.20, 10);
		// a.appendNewMatterInterval(qww, Eg - dEg, 0.25, 10);
		// a.appendNewMatterInterval(qww, Eg - dEg / 3, 0.25, 10);
		// a.appendNewMatterInterval(qwb, Eg, 0.20, 10);
		// a.appendNewMatterInterval(qww, Eg - dEg, 0.25, 10);
		// a.appendNewMatterInterval(b2w, Eg, 0.20, 10);

		ArrayList<MatterInterval> inlList = new ArrayList<MatterInterval>();

		inlList.add(new MatterInterval(2e-5, 1.4, 0.1002, diel));
		inlList.add(new MatterInterval(3e-7, 1, 0.067, diel));
		inlList.add(new MatterInterval(2e-5, 1.4, 0.1002, diel));

		a = new MatterStack(inlList);

		// double energy;
		// energy = 0.25 * Utils.ev2erg(dEg);
		// for (int i = 0; i < a.intervals.size(); i++) {
		// MatterInterval inl = a.getInterval(i);
		// System.out.println("k = " + inl.getWaveNumber(energy) + "\tWidth:" +
		// inl.getWidth() + "\txFromStart = "
		// + a.getSubStackWidth(1, i) + "\tDelta = " + a.getDelta(i, energy));
		// }
		// System.out.println(a);
		// System.out.println("Energy range array length: " +
		// a.getEnergyRange().length + " from " + a.getEnergyRange()[0] + " to "
		// + a.getEnergyRange()[a.getEnergyRange().length - 1] + " ergs");
		// System.out.println("Energy range array length: " +
		// a.getEnergyRange().length + " from " + a.getEnergyRange()[0] + " to "
		// + Utils.erg2ev(a.getEnergyRange()[a.getEnergyRange().length - 1]) +
		// " eV");
		// System.out.println("Eigen energy array length: " +
		// a.getEigenEnergy().length);
		// System.out.println("MatterStack width: " + a.getWidth());
		a.printEigenEnergy(45);
		System.out.println(a);
		System.out.println("Propagation matrix for region 2: "
				+ a.getPropagationMatrix(1, Utils.ev2erg(0.17)));
		System.out.println("Wave function coefficients for region 2:  "
				+ a.getWavefunctionCoeffs(1, Utils.ev2erg(0.17))[0] + "    "
				+ a.getWavefunctionCoeffs(2, Utils.ev2erg(0.17))[1]);
		System.out.println("Propagation matrix for region 1: "
				+ a.getPropagationMatrix(0, Utils.ev2erg(0.17)));
		System.out.println("Wave function coefficients for region 1:  "
				+ a.getWavefunctionCoeffs(0, Utils.ev2erg(0.17))[0] + "    "
				+ a.getWavefunctionCoeffs(0, Utils.ev2erg(0.17))[1]);
		a.plotPropagationCurves();
	}
}
