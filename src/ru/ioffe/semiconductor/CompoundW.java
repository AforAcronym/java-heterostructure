package ru.ioffe.semiconductor;

/**
 * 
 * @author Evgeny Shevchenko
 */
// public class CompoundW extends SemiconductorW {
public class CompoundW extends SemiconductorW
{

	// protected double a;
	// protected double c;
	// protected double a_0;
	// protected double c_0;
	// protected double a_Tch;
	// protected double c_Tch;
	// protected double a_TEC;
	// protected double c_TEC;
	// protected double energy_gap;
	// protected double energy_gap_0;
	// protected double alpha;
	// protected double beta;
	// protected double mass_elx;
	// protected double mass_elz;
	// protected double mass_hhx;
	// protected double mass_hhz;
	// protected double mass_lhx;
	// protected double mass_lhz;
	// protected double c11;
	// protected double c12;
	// protected double c13;
	// protected double c33;
	// protected double c44;
	// protected double pz13;
	// protected double pz33;
	// protected double pz15;
	// protected double polar_sp;
	// protected double permitivity;
	// protected double temperature;

	protected double	x;			// Compound
	protected double	alpha_bowing;
	protected double	beta_bowing;
	protected double	energy_gap_bowing;
	protected double	polar_sp_bowing;



	/**
	 * 
	 * @param x
	 * @param T
	 * @param A
	 * @param B
	 * @param energy_gap_bowing
	 * @param alpha_bowing
	 * @param beta_bowing
	 * @param polar_sp_bowing
	 */
	public CompoundW(double x, double T, SemiconductorW A, SemiconductorW B, double energy_gap_bowing, double alpha_bowing,
			double beta_bowing, double polar_sp_bowing)
	{

		this.a = vegard(A.a(), B.a());
		this.c = vegard(A.c(), B.c());
		this.a_0 = vegard(A.a_0(), B.a_0());
		this.c_0 = vegard(A.c_0(), B.c_0());
		this.a_Tch = vegard(A.a_Tch(), B.a_Tch());
		this.c_Tch = vegard(A.c_Tch(), B.c_Tch());
		this.a_TEC = vegard(A.a_TEC(), B.a_TEC());
		this.c_TEC = vegard(A.c_TEC(), B.c_TEC());
		this.energy_gap = vegard(A.energy_gap(), B.energy_gap(), energy_gap_bowing);
		this.energy_gap_0 = vegard(A.energy_gap_0(), B.energy_gap_0(), energy_gap_bowing);

		this.alpha_bowing = alpha_bowing;
		this.alpha = vegard(A.alpha(), B.alpha(), alpha_bowing);
		this.beta_bowing = beta_bowing;
		this.beta = vegard(A.beta(), B.beta(), beta_bowing);

		this.mass_elx = vegard(A.mass_elx(), B.mass_elx());
		this.mass_elz = vegard(A.mass_elz(), B.mass_elz());
		this.mass_hhx = vegard(A.mass_hhx(), B.mass_hhx());
		this.mass_hhz = vegard(A.mass_hhz(), B.mass_hhz());
		this.mass_lhx = vegard(A.mass_lhx(), B.mass_lhx());
		this.mass_lhz = vegard(A.mass_lhz(), B.mass_lhz());

		this.c11 = vegard(A.c11(), B.c11());
		this.c12 = vegard(A.c12(), B.c12());
		this.c13 = vegard(A.c13(), B.c13());
		this.c33 = vegard(A.c33(), B.c33());
		this.c44 = vegard(A.c44(), B.c44());

		this.pz13 = vegard(A.pz13(), B.pz13());
		this.pz33 = vegard(A.pz33(), B.pz33());
		this.pz15 = vegard(A.pz15(), B.pz15());

		this.polar_sp_bowing = polar_sp_bowing;
		this.polar_sp = vegard(A.polar_sp(), B.polar_sp(), polar_sp_bowing);

		this.permitivity = vegard(A.permitivity(), B.permitivity());
		this.temperature = T;
		this.energy_gap_bowing = energy_gap_bowing;

		this.x = x;
		this.temperature = T;
		A.setTemperature(T);
		B.setTemperature(T);
	}



	public final double vegard(double par1, double par2)
	{
		return Utils.vegard(x, par1, par2);
	}



	public final double vegard(double par1, double par2, double bow)
	{
		return Utils.vegard(x, par1, par2, bow);
	}

}
