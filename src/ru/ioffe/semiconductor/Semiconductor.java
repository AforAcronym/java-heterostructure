package ru.ioffe.semiconductor;

/**
 * @author Evgeny Shevchenko
 * @version 0.1
 */
abstract class Semiconductor {

    // Lattice info

    abstract double getLattice();


    abstract double getEnergyGap();


    abstract double getStaticPermitivity();


    abstract double[] getMassTensor();


    abstract double[] getPiezoTensor();


    abstract double[] getElasticTensor();


    abstract double getMassElA();


    abstract double getMassElB();


    abstract double getMassElC();


    abstract double getMassLHA();


    abstract double getMassLHB();


    abstract double getMassLHC();


    abstract double getMassHHA();


    abstract double getMassHHB();


    abstract double getMassHHC();

}
