/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
package ru.ioffe.semiconductor;

import org.opensourcephysics.numerics.Integral;

/**
 * Here are based some utility functions
 * 
 * @author Evgeny Shevchenko
 * @version 0.1
 * 
 */
public strictfp class Utils
{

	/**
	 * Vegard's law
	 * 
	 * @param x
	 *                compound in A(x)B(1-x)
	 * @param par1
	 *                any property of A(x)
	 * @param par2
	 *                any property of B(1-x)
	 * @param bowing
	 *                bowing factor
	 * @return the property of A(x)B(1-x)
	 */
	public static double vegard(double x, double par1, double par2, double bowing)
	{
		return x * par1 + (1 - x) * par2 - x * (1 - x) * bowing;
	}



	/**
	 * Vegard's law without parabolic member
	 * 
	 * @param x
	 *                compound in A(x)B(1-x)
	 * @param par1
	 *                any property of A(x)
	 * @param par2
	 *                any property of B(1-x)
	 * @return the property of A(x)B(1-x)
	 */
	public static double vegard(double x, double par1, double par2)
	{
		return vegard(x, par1, par2, 0);
	}



	/**
	 * Varshni equaton for energy gap dependence from temperature
	 * 
	 * @param T
	 *                temperature, K
	 * @param e0
	 *                energy gap at 0 K in eV
	 * @param aplha
	 *                eV/K
	 * @param beta
	 *                K
	 * @return energy gap, eV
	 */
	public static double varshni(double T, double e0, double alpha, double beta)
	{
		return e0 - alpha * T * T / (beta + T);
	}



	/**
	 * Debye function for lattice parameters calculation depending on
	 * material temperature
	 * 
	 * @param x
	 *                use to be temperatures ratio
	 * @param prec
	 *                precision, defines number of steps from 0 to 1
	 * @return double lattice parameter, angstrom
	 */
	public static double debye(double x, int prec)
	{
		double step = 1 / (double) prec; // Integration step (0..1)
		double[] f = new double[prec + 1]; // Function values array
		
		for (int i = 0; i < f.length; i++)
		{
			f[i] = Math.pow(i * step, 3) / (Math.exp(i * step * x) - 1);
		}
		return 3 * Integral.simpson(f, step);
	}



	/**
	 * Debye function for lattice parameters calculation depending on
	 * material temperature with default number of integration steps 10000
	 * 
	 * @param x
	 *                usually a temperatures ratio
	 * @return (double) lattice parameter, in angstroms
	 */
	public static double debye(double x)
	{
		int prec = 10000; // Default precision
		return debye(x, prec);
	}



	/**
	 * Converts electron-volts to ergs (energy in CGS)
	 * 
	 * @param e
	 *                energy in electron-volts
	 * @return
	 */
	public static double ev2erg(double e)
	{
		return e * General.EL_CHARGE_SI * 1e7;
	}



	/**
	 * Converts ergs to electron-volts
	 * 
	 * @param e
	 *                energy in ergs
	 * @return energy in electron-volts
	 */
	public static double erg2ev(double e)
	{
		return e / General.EL_CHARGE_SI * 1e-7;
	}



	/**
	 * Converts Volts per cm to CGS electric field units
	 * 
	 * @param f
	 *                electric field in V/cm
	 * @return electric field in CGS units
	 */
	public static double vpcm2CGS(double f)
	{
		return f * 1e8 / General.LIGHT_SPEED;
	}



	/**
	 * Converts CGS electric field units to Volts per cm
	 * 
	 * @param f
	 *                electric field in CGS units
	 * @return electric field in V/cm
	 */
	public static double CGS2Vpcm(double f)
	{
		return f * 1e-8 * General.LIGHT_SPEED;
	}
}
