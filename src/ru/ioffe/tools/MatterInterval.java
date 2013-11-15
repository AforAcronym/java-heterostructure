/**
 * Interval of solid material for calculating wave function in it
 */
package ru.ioffe.tools;

import ru.ioffe.semiconductor.General;
import ru.ioffe.semiconductor.Utils;

/**
 * @author Evgeny Shevchenko
 * @version 0.1
 */
public strictfp class MatterInterval extends Interval {

    private double mass;
    private double dielectricConstant;
    private double concentrationAcceptor;
    private double concentrationDonor;
    private double concentrationElectrons;
    private double concentrationHoles;


    /**
     * @param x_start
     * @param x_end
     * @param y
     */
    public MatterInterval(double x_start, double x_end, double y) {
        super(x_start, x_end, y);
        mass = 1;
        dielectricConstant = 1;
    }


    /**
     * @param width
     * @param y
     * @param mass
     * @param diel
     */
    public MatterInterval(double width, double y, double mass, double diel) {
        super(0, width, y);
        this.mass = mass;
        this.dielectricConstant = diel;
    }


    /**
     * @param xstart
     * @param xend
     * @param y
     * @param mass
     * @param diel
     */
    public MatterInterval(double xstart, double xend, double y, double mass, double diel) {
        super(xstart, xend, y);
        this.mass = mass;
        this.dielectricConstant = diel;
    }


    @Override
    public String toString() {
        return String.format("[%.2e, %.2e]\tY: %.3f\tmass: %.4f\tdiel: %.1f\t%s", x_start, x_end, Utils.erg2ev(y), mass,
                dielectricConstant, label);
    }


    public double getAcceptorConcentration() {
        return concentrationAcceptor;
    }


    public void setAcceptorConcentration(double acceptorConcentration) {
        this.concentrationAcceptor = acceptorConcentration;
    }


    public double getDonorConcentration() {
        return concentrationDonor;
    }


    public void setDonorConcentration(double donorConcentration) {
        this.concentrationDonor = donorConcentration;
    }


    public double getElectronsConcentration() {
        return concentrationElectrons;
    }


    public void setElectronsConcentration(double electronsConcentration) {
        this.concentrationElectrons = electronsConcentration;
    }


    public double getHolesConcentration() {
        return concentrationHoles;
    }


    public void setHolesConcentration(double holesConcentration) {
        this.concentrationHoles = holesConcentration;
    }


    public double getMass() {
        return mass;
    }


    public double getDiel() {
        return dielectricConstant;
    }


    /**
     * Get wave number for the interval and chosen energy
     *
     * @param energy
     * @return
     */
    public Complex getWaveNumber(double energy) {
        return Complex.pow(energy - y, 0.5).times(Math.sqrt(2 * mass * General.EL_MASS) / General.PLANCK_HBAR_ERGS);
    }


    /**
     * Main method bleat
     *
     * @param args
     */
    public static void main(String[] args) {

        double x0 = 0;
        double step = 0.5;
        double v = 3;
        double m = 0.1;
        double d = 10;
        double energy = v / 2;

        MatterInterval i1 = new MatterInterval(x0, x0 + step, v, m, d);
        MatterInterval i2 = new MatterInterval(i1.getXEnd(), i1.getXEnd() + step, 0 / 2, m * 0.9, d * 1.1);
        MatterInterval i3 = new MatterInterval(i2.getXEnd(), i2.getXEnd() + step, v, m, d);

        System.out.println(i1);
        System.out.println(i2);
        System.out.println(i3);
        System.out.println("k1 = " + i1.getWaveNumber(energy) + " cm^-1");
        System.out.println("k2 = " + i2.getWaveNumber(energy) + " cm^-1");
        System.out.println("k3 = " + i3.getWaveNumber(energy) + " cm^-1");

        // System.out.println(1.0 / 0); // Infinity
        // (Double.POSITIVE_INFINITY)

    }

}
