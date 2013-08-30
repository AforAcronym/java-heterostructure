package ru.ioffe.semiconductor;

/**
 * 
 * @author Evgeny Shevchenko
 * @version 0.1
 */
abstract class Semiconductor {

	// Lattice info
	abstract double getA();

	abstract double getB();

	abstract double getC();

	abstract double getCellVolume();

	abstract double getEnergyGap();

	abstract double getStaticPermitivity();

	abstract double[] getMassTensor();

	abstract double[] getPiezoTensor();

	abstract double[] getElasticTensor();

	abstract double get_mass_el_a();

	abstract double get_mass_el_b();

	abstract double get_mass_el_c();

	abstract double get_mass_lh_a();

	abstract double get_mass_lh_b();

	abstract double get_mass_lh_c();

	abstract double get_mass_hh_a();

	abstract double get_mass_hh_b();

	abstract double get_mass_hh_c();

}
