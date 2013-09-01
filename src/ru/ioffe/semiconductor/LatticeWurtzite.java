/**
 * Wurtzite lattice class
 */
package ru.ioffe.semiconductor;

/**
 * @author Evgeny Shevchenko
 * 
 */
public strictfp class LatticeWurtzite extends Lattice {

	// Lattice parameters
	private double a, c;

	public LatticeWurtzite(double a, double c) {
		this.a = a;
		this.c = c;
	}
    
    protected double a_0, c_0;  // Lattice constants at T = 0
    protected double a_TEC, c_TEC; // Temperature expansion coefficients, 1/K
    protected double a_Tch, c_Tch; // Characteristic (Debye) temperature, K
    
  
	/**
	 * Cell volume
	 * 
	 * @return volume in nm^3
	 */
	public double getCellVolume() {
		// Math.pow(3., 0.5) = 1.73205080757
		return 0.5 * 3 * 1.73205080757 * a * a * c;
	}


	@Override
	public void setTemperature(double T) {
		this.temperature = T;
        a = latticePar(T, a_0, a_TEC, a_Tch);
        c = latticePar(T, c_0, c_TEC, c_Tch);
	}

}
