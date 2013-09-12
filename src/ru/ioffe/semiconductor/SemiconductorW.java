package ru.ioffe.semiconductor;

/**
 * Semiconductor with wurtzite lattice type
 * 
 * @author Evgeny Shevchenko
 */
public strictfp class SemiconductorW
{

	protected double	a, c;	// Lattice constants, angstrom
	protected double	a_0, c_0;	// Lattice constants at T = 0
	protected double	a_TEC, c_TEC;	// Temperature expansion coefficients
	protected double	a_Tch, c_Tch;	// Characteristic (Debye) temperature
	protected double	energy_gap, energy_gap_0;	// Energy gap at 300 and 0 K
								// correspondingly
	protected double	alpha, beta;			// Coefficients for Varshni equation
	protected double	mass_elx, mass_elz;		// Effective masses of electrons
	protected double	mass_hhx, mass_hhz;		// Effective masses of heavy holes
	protected double	mass_lhx, mass_lhz;		// Effective masses of light holes
	protected double	c11, c12, c13, c33, c44;	// Elastic coeffients
	protected double	pz13, pz33, pz15;		// Piezoelectric tensor components
	protected double	polar_sp;			// Spontaneous polarization
	protected double	permitivity;			// Static dielectric constant
	protected double	temperature;			// Temperature



	public SemiconductorW()
	{
		// Just an empty constructor needed for inheritance
	}



	/**
	 * 
	 * @param a
	 *                Lattice constant a, angstrom
	 * @param c
	 *                Lattice constant c, angstrom
	 * @param a_0
	 *                Lattice constant a at T=0, angstroms
	 * @param c_0
	 *                Lattice constant c at T=0, angstroms
	 * @param a_Tch
	 *                Characteristic (Debye) temperature for parameter a, K
	 * @param c_Tch
	 *                Characteristic (Debye) temperature for parameter c, K
	 * @param a_TEC
	 *                Temperature expansion coefficients for parameter a, 1/K
	 * @param c_TEC
	 *                Temperature expansion coefficients for parameter c, 1/K
	 * @param energy_gap
	 *                Energy gap, eV
	 * @param energy_gap_0
	 *                Eergy gap at T=0, eV
	 * @param alpha
	 *                Varshni equation coefficient, eV/K
	 * @param beta
	 *                Varshni equation coefficient, K
	 * @param mass_elx
	 *                Electron effective mass in x direction, in units of free electron mass
	 * @param mass_elz
	 *                Electron effective mass in z direction, in units of free electron mass
	 * @param mass_hhx
	 *                Heavy hole effective mass in x direction, in units of free electron mass
	 * @param mass_hhz
	 *                Heavy hole effective mass in z direction, in units of free electron mass
	 * @param mass_lhx
	 *                Light hole effective mass in x direction, in units of free electron mass
	 * @param mass_lhz
	 *                Light hole effective mass in z direction, in units of free electron mass
	 * @param c11
	 *                Elastic coeffient, GPa
	 * @param c12
	 *                Elastic coeffient, GPa
	 * @param c13
	 *                Elastic coeffient, GPa
	 * @param c33
	 *                Elastic coeffient, GPa
	 * @param c44
	 *                Elastic coeffient, GPa
	 * @param pz13
	 *                Piezoelectric tensor components, C/cm^2
	 * @param pz33
	 *                Piezoelectric tensor components, C/cm^2
	 * @param pz15
	 *                Piezoelectric tensor components, C/cm^2
	 * @param polar_sp
	 *                Spontaneous polarization, C/cm^2
	 * @param permitivity
	 *                Dielectric constant
	 * @param T
	 *                Temperature, K
	 */
	public SemiconductorW(double a, double c, double a_0, double c_0, double a_Tch, double c_Tch, double a_TEC, double c_TEC,
			double energy_gap, double energy_gap_0, double alpha, double beta, double mass_elx, double mass_elz,
			double mass_hhx, double mass_hhz, double mass_lhx, double mass_lhz, double c11, double c12, double c13,
			double c33, double c44, double pz13, double pz33, double pz15, double polar_sp, double diel, double T)
	{
		this.a = a;
		this.c = c;
		this.a_0 = a_0;
		this.c_0 = c_0;
		this.a_Tch = a_Tch;
		this.c_Tch = c_Tch;
		this.a_TEC = a_TEC;
		this.c_TEC = c_TEC;
		this.energy_gap = energy_gap;
		this.energy_gap_0 = energy_gap_0;
		this.alpha = alpha;
		this.beta = beta;
		this.mass_elx = mass_elx;
		this.mass_elz = mass_elz;
		this.mass_hhx = mass_hhx;
		this.mass_hhz = mass_hhz;
		this.mass_lhx = mass_lhx;
		this.mass_lhz = mass_lhz;
		this.c11 = c11;
		this.c12 = c12;
		this.c13 = c13;
		this.c33 = c33;
		this.c44 = c44;
		this.pz13 = pz13;
		this.pz33 = pz33;
		this.pz15 = pz15;
		this.polar_sp = polar_sp;
		this.permitivity = diel;
		this.temperature = T;
		setTemperature(T);

	}



	/**
	 * This method sets temperature dependent values of bulk semiconductor
	 * 
	 * @param T
	 *                temperature
	 */
	public final void setTemperature(double T)
	{
		this.temperature = T;
		energy_gap = Utils.varshni(T, energy_gap_0, alpha, beta);
		a = latticePar(T, a_0, a_TEC, a_Tch);
		c = latticePar(T, c_0, c_TEC, c_Tch);
	}



	/**
	 * Calculation of lattice parameter
	 * 
	 * @param T
	 *                The temperature
	 * @param par0
	 *                Lattice parameter at T=0
	 * @param TEC
	 *                Thermal expansion coefficient
	 * @param Tch
	 *                Debye temperature
	 * @return Lattice parameter, angstroms
	 */
	private double latticePar(double T, double par0, double TEC, double Tch)
	{
		if (temperature == 0)
		{
			return par0;
		}
		else
		{
			return par0 * (1 + TEC * Tch * Utils.debye(Tch / T));
		}
	}



	/**
	 * Unit cell volume
	 * 
	 * @return angstrom^3
	 */
	public double cellVolume()
	{
		return 0.5 * Math.pow(3., 1.5) * a * a * c;
	}



	public double reducedMass(char direction)
	{
		if (direction == 'x' || direction == 'a')
		{
			return 1 / (1 / mass_elx + 1 / mass_hhx);
		}
		else
		{
			return 1 / (1 / mass_elz + 1 / mass_hhz);
		}
	}



	/**
	 * Exciton Bohr radius
	 * 
	 * @param direction
	 *                'x' or 'a' in order not along 'c' or 'z'
	 * @return
	 */
	public double excBohrRadius(char direction)
	{
		double rmass = reducedMass(direction);
		return permitivity * Math.pow(General.PLANCK_HBAR_ERGS, 2)
				/ (rmass * General.EL_MASS * Math.pow(General.EL_CHARGE_CGS, 2));

	}



	/**
	 * Exciton binding energy
	 * 
	 * @param direction
	 *                'x' or 'a' in order not along 'c' or 'z'
	 * @param n
	 *                number of exciton energy level
	 * @return double
	 */
	public double excBindEnergy(char direction, int n)
	{
		return Math.pow(General.EL_CHARGE_CGS, 2) / 2 / permitivity / excBohrRadius(direction) / n / n * 1e-7
				/ General.EL_CHARGE_SI;
	}



	public double a()
	{
		return a;
	}



	public double c()
	{
		return c;
	}



	public double a_0()
	{
		return a_0;
	}



	public double c_0()
	{
		return c_0;
	}



	public double a_TEC()
	{
		return a_TEC;
	}



	public double c_TEC()
	{
		return c_TEC;
	}



	public double a_Tch()
	{
		return a_Tch;
	}



	public double c_Tch()
	{
		return c_Tch;
	}



	public double energy_gap()
	{
		return energy_gap;
	}



	public double energy_gap_0()
	{
		return energy_gap_0;
	}



	public double alpha()
	{
		return alpha;
	}



	public double beta()
	{
		return beta;
	}



	public double mass_elx()
	{
		return mass_elx;
	}



	public double mass_elz()
	{
		return mass_elz;
	}



	public double mass_hhx()
	{
		return mass_hhx;
	}



	public double mass_hhz()
	{
		return mass_hhz;
	}



	public double mass_lhx()
	{
		return mass_lhx;
	}



	public double mass_lhz()
	{
		return mass_lhz;
	}



	public double c11()
	{
		return c11;
	}



	public double c12()
	{
		return c12;
	}



	public double c13()
	{
		return c13;
	}



	public double c33()
	{
		return c33;
	}



	public double c44()
	{
		return c44;
	}



	public double pz13()
	{
		return pz13;
	}



	public double pz33()
	{
		return pz33;
	}



	public double pz15()
	{
		return pz15;
	}



	public double polar_sp()
	{
		return polar_sp;
	}



	public double permitivity()
	{
		return permitivity;
	}



	public double temperature()
	{
		return temperature;
	}
}
