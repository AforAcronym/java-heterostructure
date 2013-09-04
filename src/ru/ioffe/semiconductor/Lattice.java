package ru.ioffe.semiconductor;

abstract class Lattice {

	// Lattice parameters at desired temperature T
	// protected double a, b, c;
	protected double a;

	// Lattice parameters at T=0
	// protected double a0, b0, c0;
	protected double a0;

	// Temperature in Kelvins
	protected double temperature;

	// TEC - Temperature expansion coefficients for parameter a, 1/K

	protected double a_TEC;
	// protected double a_TEC, b_TEC, c_TEC;

	// Tch - Characteristic (Debye) temperature for parameter c, K
	// protected double a_Tch, b_Tch, c_Tch;
	protected double a_Tch;

	/**
	 * This method sets temperature dependent values of bulk semiconductor
	 * 
	 * @param T
	 *            temperature in Kelvins
	 */
	public abstract void setTemperature(double T);

	/**
	 * Calculation of lattice parameter
	 * 
	 * @param T
	 *            The temperature
	 * @param par_0
	 *            Lattice parameter at T = 0 K
	 * @param TEC
	 *            Thermal expansion coefficient
	 * @param Tch
	 *            Debye temperature
	 * @return Lattice parameter
	 */
	protected double latticePar(double T, double par_0, double TEC, double Tch) {
		if (T == 0) {
			return par_0;
		} else {
			return par_0 * (1 + TEC * Tch * Utils.debye(Tch / T));
		}
	}

	public double getA() {
		return a;
	}
//
//	public double getB() {
//		return b;
//	}
//
//	public double getC() {
//		return c;
//	}

	public double getCellVolume() {
		return a * a * a;
	}
}
