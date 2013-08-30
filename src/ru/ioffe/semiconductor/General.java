package ru.ioffe.semiconductor;

/**
 * General matter properties
 * @author Evgeny Shevchenko
 *
 */
public abstract  strictfp class General {

    /** 
     * Vacuum permitivity (electric constant) in Farads per meter
     */
    public static final double VAC_PERMIT = 8.854187817620e-12; 
    
    /**
     * Light speed in vacuum, cm / sec
     */
    public static final double LIGHT_SPEED = 29979245800.;	
    
    /** 
     * Electron charge in SGC, L^1.5 * M^0.5 / T 
     */
    public static final double EL_CHARGE_CGS = 4.80325e-10; 
    
    /** 
     * Electron charge in SI, Coloumb
     */
    public static final double EL_CHARGE_SI = 1.602176565e-19; 
    
    /** 
     * Free electron mass, g
     */
    public static final double EL_MASS = 9.10938188e-28; 
    
    /**
     * Planck constant, J*s
     */
    public static final double PLANCK_H_JS = 6.62606957e-34; 
    
    /** 
     * Planck constant, erg*s
     */
    public static final double PLANCK_H_ERGS = PLANCK_H_JS * 1e7; 
    
    /** 
     * Planck constant, eV*s
     */
    public static final double PLANCK_H_EVS = PLANCK_H_JS / EL_CHARGE_CGS; 
    
    /** 
     * Reduced planck constant, J*s
     */
    public static final double PLANCK_HBAR_JS = PLANCK_H_JS / 2 / Math.PI; 
    
    /** 
     * Reduced planck constant, erg*s
     */
    public static final double PLANCK_HBAR_ERGS = PLANCK_H_ERGS / 2 / Math.PI; 
    
    /** 
     * Reduced planck constant, eV*s
     */
    public static final double PLANCK_HBAR_EVS = PLANCK_H_EVS / 2 / Math.PI; 
    
    /** 
     * Rydberg, eV
     */
    public static final double RY_EV = 13.60569253; 
}
