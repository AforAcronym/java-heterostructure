/**
 *
 */
package ru.ioffe.tools;

import org.opensourcephysics.frames.PlotFrame;
import ru.ioffe.semiconductor.Utils;

import javax.swing.*;
import java.util.ArrayList;

/**
 * @author Evgeny Shevchenko
 */
public strictfp class MatterStackGhatak {

    // 0.0001 eV in ergs
    private double ENERGY_STEP = Utils.ev2erg(0.001);

    // point per 1 nm
    private final int X_POINTS_DENSITY = 40;

    private double[] eigenEnergy = null;
    private double[] energyRange = null;
    private ArrayList<MatterInterval> intervals = new ArrayList<MatterInterval>();
    private double y_shift = 0;
    private boolean fictitiousRegionReady = false;

    private boolean VERBOSE = false;


    // private boolean VERBOSE = true;

    /**
     * *************************************************************************************************************************
     */
    public MatterStackGhatak() {
    }


    /**
     * *************************************************************************************************************************
     */
    public MatterStackGhatak(ArrayList<MatterInterval> intervals) {
        this.intervals = intervals;
        addFictitiousInterval();
        process();
    }


    /**
     * *************************************************************************************************************************
     */
    @Override
    public String toString() {
        StringBuffer s = new StringBuffer();
        int counter = 0;
        for (MatterInterval mi : intervals) {
            s.append("" + (counter++) + " : ");
            s.append(mi.toString() + "\n");
        }
        return s.toString();
    }


    /**
     * *************************************************************************************************************************
     */
    public MatterInterval getInterval(int index) {
        return intervals.get(index);
    }


    /**
     * *************************************************************************************************************************
     */
    private void addFictitiousInterval() {
        if (intervals.size() == 0) {
            throw new IndexOutOfBoundsException("There must be at least one interval to add fictitious one");
        }
        // TODO take mass of the QWs for the fictitious layer
        // intervals.add(0, new MatterInterval(0, 0, 0,
        // intervals.get(0).getMass(), intervals.get(0)
        // .getDiel()));
        intervals.add(0, new MatterInterval(0, 0, 0, intervals.get(0).getMass(), intervals.get(0).getDiel()));
        // intervals.add(0, new
        // MatterInterval(-intervals.get(0).getWidth(), 0, 0,
        // intervals.get(0).getMass(), intervals.get(0).getDiel()));
        getFictitiousInterval().setLabel("Fictitious");
        fictitiousRegionReady = true;
    }


    /**
     * *************************************************************************************************************************
     * Access to fictitious interval
     *
     * @return MatterInterval instance of the fictitious interval that is always 0th in the MatterStack
     */
    private MatterInterval getFictitiousInterval() {
        return intervals.get(0);
    }


    /**
     * *************************************************************************************************************************
     * Add an interval to the end of the stack
     *
     * @param element MatterInterval instance
     */
    public void append(MatterInterval element) {
        add(intervals.size(), element);
    }


    /**
     * *************************************************************************************************************************
     * Add an interval to the chosen position of the stack. Index should start from 1 in order not to affect fictitious 0th
     * interval
     *
     * @param index   position
     * @param element MatterInterval instance
     */
    public void add(int index, MatterInterval element) {
        if (index < 1 && fictitiousRegionReady) {
            throw new IndexOutOfBoundsException("Index must be from 1 to " + intervals.size() + " to add to MatterStack, but "
                    + index + " inserted");
        }
        processBack();
        intervals.add(index, element);
        if (!fictitiousRegionReady) {
            addFictitiousInterval();
        }
        process();
    }


    /**
     * *************************************************************************************************************************
     * Add a new MatterInterval region, specifying its essential properties for constructor
     *
     * @param width width of the interval
     * @param y     potential value in electron-volts
     * @param mass  mass of a particle in units of free electron mass
     * @param diel  dielectric constant of the interval
     */
    public void appendNewMatterInterval(double width, double y, double mass, double diel) {
        double xstart = getWidth();
        add(intervals.size(), new MatterInterval(xstart, xstart + width, Utils.ev2erg(y), mass, diel));
    }


    /**
     * *************************************************************************************************************************
     * Process potential to reduce the lowest level to zero level
     */
    private void process() {
        if (intervals.size() < 2) {
            throw new IndexOutOfBoundsException("Processing requires at least two intervals");
        }
        // addressing to the first non-fictitious interval
        double y_min = intervals.get(1).getY();
        MatterInterval inl;

        for (int i = 1; i < intervals.size(); i++) {
            inl = intervals.get(i);
            if (inl.getY() < y_min) {
                y_min = inl.getY();
            }
        }
        y_shift = y_min;

        for (int i = 1; i < intervals.size(); i++) {
            inl = intervals.get(i);
            inl.setY(inl.getY() - y_shift);
        }

        getFictitiousInterval().setY(0);
        // System.out.println("After process:\n" + this);
    }


    /**
     * *************************************************************************************************************************
     * Returns initial potential of intervals for correct processing
     */
    private void processBack() {
        if (intervals == null || intervals.size() == 0) {
            return;
        }
        getFictitiousInterval().setY(y_shift);
        for (int i = 1; i < intervals.size(); i++) {
            getInterval(i).setY(getInterval(i).getY() + y_shift);
        }
        // System.out.println("After processBack:\n" + this);
    }


    /**
     * *************************************************************************************************************************
     */
    public Complex getDelta(int index, double energy) {
        MatterInterval inl = intervals.get(index);
        if (index == 0) {
            // Fictitious interval
            return inl.getWaveNumber(energy).times(-inl.getWidth());
        } else if (index == 1) {
            return new Complex(0.0, 0.0);
        } else {
            // FIXME Check sum of the intervals widths from 0 to
            // index.
            // The intervals must not have equal widths, for now the
            // opposite is
            // assumed -- TEST THE LINE BELOW
            return inl.getWaveNumber(energy).times(getSubStackWidth(1, index));
        }
    }


    /**
     * *************************************************************************************************************************
     * Takes index of preceding interval (between index and index+1) and energy value
     *
     * @return ComplexMatrix
     * @throws Exception
     */
    public ComplexMatrix getLocalPropagationMatrix(int index, double energy) {

        MatterInterval inl = intervals.get(index);
        MatterInterval next = intervals.get(index + 1);

        Complex k_ratio;
        // Double.POSITIVE_INFINITY is handled here and it assigns to real part of Complex if appears
        k_ratio = Complex.sqrt((energy - next.getY()) / (energy - inl.getY()) * inl.getMass() / next.getMass());

        Complex kcp = k_ratio.add(1); // 1 + km_ratio
        Complex kcm = k_ratio.reverseSign().add(1); // 1 - km_ratio
        Complex kd = inl.getWaveNumber(energy).times(inl.getWidth());

        if (VERBOSE) {
            System.out.println("k_ratio = " + k_ratio);
            System.out.println("kcp = " + kcp);
            System.out.println("kcm = " + kcm);
            System.out.println("kd = " + kd);
        }

        Complex i = new Complex(0, 1);
        Complex ni = i.reverseSign();
        try {
            return new ComplexMatrix(new Complex[][]{
                    {kcp.times(Complex.exp(ni.times(kd))), kcm.times(Complex.exp(ni.times(kd)))},
                    {kcm.times(Complex.exp(i.times(kd))), kcp.times(Complex.exp(i.times(kd)))}}).multiplyEach(0.5);
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }


    /**
     * *************************************************************************************************************************
     * Get array of complex coefficients of exponential 1-dimensional wave function
     *
     * @param index
     * @param energy
     * @return
     */
    public ComplexMatrix getPropagationMatrix(int index, double energy) {
        ComplexMatrix mx = ComplexMatrix.identityMatrix(2);
        for (int i = index; i < intervals.size() - 1; i++) {
            try {
                mx = mx.multiply(getLocalPropagationMatrix(i, energy));
                if (VERBOSE) {
                    System.out.println("> > > In getIntervalWavefunctionCoeffs(...)");
                    System.out.println("Matrix on index " + index + ", energy " + energy + ", interation " + i);
                    System.out.println(mx);
                    System.out.println("Current propagation matrix = " + getLocalPropagationMatrix(i, energy));
                }
            } catch (NullPointerException e) {
                System.out.println("Null received in ComplexMatrix multiplication in getIntervalWavefunctionCoeffs(...)");
                e.printStackTrace();
            } catch (Exception e) {
                System.out.println("Something wrong on ComplexMatrix multiplication in getIntervalWavefunctionCoeffs(...)");
                e.printStackTrace();
            }
        }
        return mx;
    }


    /**
     * *************************************************************************************************************************
     * Returns array of energy range from zero potential to the second maximum so that at least one QW is taken into account
     */
    public double[] getEnergyRange() {
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
        int rangesize = (int) (max2 / ENERGY_STEP);
        if (VERBOSE) {
            System.out.println("Energy max:  " + Utils.erg2ev(max));
            System.out.println("Energy max2: " + Utils.erg2ev(max2));
            System.out.println("Energy range size: " + rangesize);
        }
        energyRange = new double[rangesize];
        for (int i = 0; i < rangesize; i++) {
            energyRange[i] = i * ENERGY_STEP;
        }
        return energyRange;
    }


    /**
     * *************************************************************************************************************************
     * Value of propagation curve. Propagation maximums correspond to eigen states
     *
     * @param index  index of a region/interval
     * @param energy energy in ergs
     * @return
     */
    public double getPropagationValue(int index, double energy) {
        if (VERBOSE) {
            System.out.println("Propagation value for index" + index + " and energy " + energy);
            // System.out.println("" + getPropagationMatrix(index,
            // energy)[0].abs());
            // System.out.println("" + getPropagationMatrix(0,
            // energy)[0].abs());
        }
        return Math.pow(getWavefunctionCoeffs(index, energy)[0].abs() / getWavefunctionCoeffs(0, energy)[0].abs(), 2);
    }


    /**
     * *************************************************************************************************************************
     * Derivative of the propagation curve
     *
     * @param index  index of a region/interval
     * @param energy energy in ergs
     * @return
     */
    public double getPropagationValueDerivative(int index, double energy) {
        return (getPropagationValue(index, energy + Function.CALC_TOLERANCE) - getPropagationValue(index, energy
                - Function.CALC_TOLERANCE))
                / 2 / Function.CALC_TOLERANCE;
    }



	/*
     * Needed only for Newton-Raphson method public double getPropagationValueDerivative2(int index, double energy) { return
	 * (getPropagationValueDerivative(index, energy + Function.CALC_TOLERANCE) - getPropagationValueDerivative( index, energy -
	 * Function.CALC_TOLERANCE)) / 2 / Function.CALC_TOLERANCE; }
	 */

    /**
     * *************************************************************************************************************************
     * Log of the value of propagation curve. Propagation maximums correspond to eigen states
     *
     * @param index  index of a region/interval
     * @param energy energy in ergs
     * @return
     */
    private double getPropagationValueLog(int index, double energy) {
        return Math.log(getPropagationValue(index, energy));
    }


    /**
     * *************************************************************************************************************************
     * Derivative of the Log of propagation curve
     *
     * @param index  index of a region/interval
     * @param energy energy in ergs
     * @return
     */
    @SuppressWarnings("unused")
    private double getPropagationValueLogDerivative(int index, double energy) {
        return (getPropagationValueLog(index, energy + Function.CALC_TOLERANCE) - getPropagationValueLog(index, energy
                - Function.CALC_TOLERANCE))
                / 2 / Function.CALC_TOLERANCE;
    }


    /**
     * *************************************************************************************************************************
     * Get propagation energy in an interval between two energy values by bisection method
     *
     * @param index      index of region/interval
     * @param energy_min minimal energy in ergs for the range to search in
     * @param energy_max maximum energy in ergs for the range to search in
     * @return energy value in ergs if maximum is found
     */
    public double getPropagationEnergy(int index, double energy_min, double energy_max) {
        System.out.println("ENTERING BISECTION...");
        double energy = (energy_min + energy_max) / 2;
        // double dprop = getPropagationValueLogDerivative(index, energy);
        double dprop = getPropagationValueDerivative(index, energy);

        while (Math.abs(dprop) > Function.CALC_TOLERANCE) {
            System.out.println("BSF: " + dprop + "\tenergy_min: " + energy_min + "\tenergy: " + energy + "\tenergy_max: "
                    + energy_max);
            if (Math.abs(energy / (energy_min + energy_max) - 0.5) < Function.CALC_TOLERANCE * 0.001) // XXX Magic
            {
                return energy;
            }

            // if (dprop * getPropagationValueLogDerivative(index, energy_min) > 0)
            if (dprop * getPropagationValueDerivative(index, energy_min) > 0) {
                energy_min = energy;
            }
            // else if (dprop * getPropagationValueLogDerivative(index, energy_max) > 0)
            else if (dprop * getPropagationValueDerivative(index, energy_max) > 0) {
                energy_max = energy;
            }

            System.out.println("Energy difference = " + (energy - (energy_min + energy_max) / 2));
            energy = (energy_min + energy_max) / 2;
            // dprop = getPropagationValueLogDerivative(index, energy);
            dprop = getPropagationValueDerivative(index, energy);
        }
        return energy;
    }


    /**
     * *************************************************************************************************************************
     * Get eigen energies for an interval
     *
     * @param index of an interval
     * @return double[] of eigen energies in ergs
     */
    public double[] findIntervalEigenEnergies(int index) {
        ArrayList<Double> energyList = new ArrayList<Double>();
        double[] eRange = getEnergyRange();
        for (int i = 1; i < eRange.length - 1; i++) {
            // if (VERBOSE) {
            // System.out.println("Index of interval: " + index + "    energy: " + eRange[i] + " erg, " +
            // Utils.erg2ev(eRange[i]) + " eV");
            // System.out.println("PropagationValueDerivative(" +
            // (index - 1) + ", " + eRange[i] + ") = " + getPropagationValueDerivative(index - 1,
            // eRange[i]));
            // System.out.println("PropagationValueDerivative(" +
            // (index) + ", " + eRange[i] + ") = " + getPropagationValueDerivative(index, eRange[i]));
            // System.out.println("PropagationValueDerivative(" + (index + 1) + ", " + eRange[i] + ") = "
            // + getPropagationValueDerivative(index + 1, eRange[i]));
            // }

            if (Double.isNaN(getPropagationValueDerivative(index, eRange[i - 1]))
                    || Double.isNaN(getPropagationValueDerivative(index, eRange[i]))
                    || Double.isNaN(getPropagationValueDerivative(index, eRange[i + 1]))) {
                System.out.println("NaNs detected!");
                continue;
            } else if (getPropagationValue(index, eRange[i]) > getPropagationValue(index, eRange[i - 1])
                    && getPropagationValue(index, eRange[i]) > getPropagationValue(index, eRange[i + 1])
                    && getPropagationValueDerivative(index, eRange[i - 1]) * getPropagationValueDerivative(index, eRange[i + 1]) < 0) {
                // energyList.add(eRange[i]);
                energyList.add(getPropagationEnergy(index, eRange[i - 1], eRange[i + 1]));
                // System.out.println("FOUND " + eRange[i] + "\tCALCULATED " + getPropagationEnergy(index,
                // eRange[i - 1], eRange[i + 1]));
            }
        }
        double[] energyArray = new double[energyList.size()];
        for (int i = 0; i < energyList.size(); i++) {
            energyArray[i] = Double.valueOf(energyList.get(i));
        }
        return energyArray;
    }


    /**
     * *************************************************************************************************************************
     * Build a plot of quantum wells propagation curves
     */
    public void plotPropagationCurves() {
        double[] eRange = getEnergyRange();
        PlotFrame frame = new PlotFrame("Energy", "Propagation", "QW Propagation Curves");
        for (int i = 1; i < intervals.size() - 1; i++) {
            if (intervals.get(i).getY() < intervals.get(i - 1).getY() && intervals.get(i).getY() < intervals.get(i + 1).getY()) {
                for (double e : eRange) {
                    frame.append(i, Utils.erg2ev(e), getPropagationValueLog(i, e));
                }
                frame.setMarkerShape(i, 0);
            }
        }
        frame.setSize(600, 600);
        frame.setConnected(true);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }


    /**
     * *************************************************************************************************************************
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

        // Among all intervals we are looking for those which are truly quantum wells, they bring all the information.
        for (int i = 1; i < intervals.size() - 1; i++) {
            if (intervals.get(i).getY() <= intervals.get(i - 1).getY() && intervals.get(i).getY() < intervals.get(i + 1).getY()) {
                for (double e : findIntervalEigenEnergies(i)) {
                    if (!eigenEnergyList.contains(Double.valueOf(e))) {
                        eigenEnergyList.add(new Double(e));
                    }
                }
            }
        }
        eigenEnergy = new double[eigenEnergyList.size()];
        for (int i = 0; i < eigenEnergy.length; i++) {
            eigenEnergy[i] = Double.valueOf(eigenEnergyList.get(i));
        }
        return eigenEnergy;
    }


    /**
     * *************************************************************************************************************************
     * Print a column of eigen energy values of the stack
     *
     * @param number number of states starting with the lowest one (not implemented yet)
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
     * *************************************************************************************************************************
     * Get the whole MatterStack width starting from index=1 (thus excluding the fictitious interval) to the last one
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
     * *************************************************************************************************************************
     * Get width of specified interval. Note that indexes start with 1 due to fictitious interval (0th)
     *
     * @param index
     * @return the interval width
     */
    public double getIntervalWidth(int index) {
        return intervals.get(index).getWidth();
    }


    /**
     * *************************************************************************************************************************
     * Get width of a part of the stack including both intervals
     *
     * @param index1 number of starting interval
     * @param index2 number of ending interval (included)
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
     * *************************************************************************************************************************
     * Get wave function coefficients for chosen region and energy
     *
     * @param index
     * @param energy
     * @return
     */
    public Complex[] getWavefunctionCoeffs(int index, double energy) {
        ComplexMatrix mx = getPropagationMatrix(index, energy);
        if (index == 1) {
            return new Complex[]{Complex.fromDouble(0), mx.getElem(1, 0)};
        }
        return new Complex[]{mx.getElem(0, 0), mx.getElem(1, 0)};
    }


    /**
     * *************************************************************************************************************************
     * Get wave function for a state
     *
     * @param energy
     * @return double[] array of a wave function values
     */
    public double[][] getWaveFunction(double energy) {

        // Number of dots in the whole MatterStack
        int numberOfDots = (int) (getWidth() * X_POINTS_DENSITY * 1000000000L);

        // TODO try-catch or think on long-or-int problem
        double[] waveFunction = new double[numberOfDots];
        double[] xAxis = new double[numberOfDots];
        double[] xLocal;

        int counter = 0;

        Complex[] coef;
        Complex jp = new Complex(0, 1);  // i
        Complex jm = new Complex(0, -1); // -i

        for (int k = 0; k < intervals.size(); k++) {
            MatterInterval inl = intervals.get(k);
            xLocal = inl.getGrid(X_POINTS_DENSITY * 1000000000L); // XXX Magic
            for (int i = 0; i < xLocal.length; i++) {

                coef = getWavefunctionCoeffs(k, energy);
                waveFunction[i + counter] = (coef[0].times(Complex.exp(jp.times(inl.getWaveNumber(energy).times(xLocal[i])
                        .subtract(getDelta(k, energy)))))).add(
                        coef[1].times(Complex.exp(jm.times(inl.getWaveNumber(energy).times(xLocal[i])
                                .subtract(getDelta(k, energy)))))).abs();

                xAxis[i + counter] = xLocal[i];

            }
            counter += xLocal.length;
        }
        return new double[][]{xAxis, waveFunction};
    }


    /**
     * *************************************************************************************************************************
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
        return new double[][]{x, y};
    }


    // private double[] reduceWaveFunction() {}

    /**
     * *************************************************************************************************************************
     *
     * @param args
     */
    @SuppressWarnings("unused")
    public static void main(String[] args) {

        double Eg = 4.0; // eV
        double dEg = 0.1; // eV
        double qww = 2e-7; // cm
        double qwb = qww / 2; // cm
        double b1w = qww * 1e2; // cm
        double b2w = b1w; // cm
        double diel = 10; // no units
        MatterStackGhatak a = new MatterStackGhatak();

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

        // appendNewMatterInterval(double width, double y, double mass, double diel)
        // width in cm, 1e-7 cm = 1 nm
        // y is energy in eV
        // mass is in free electron mass units
        // dielectric constant without units
        a.appendNewMatterInterval(200 * 1e-7, 1.4, 0.1002, diel);
        a.appendNewMatterInterval(3e-7, 1, 0.067, diel);
        a.appendNewMatterInterval(200 * 1e-7, 1.4, 0.1002, diel);

        // double energy;
        // energy = 0.25 * Utils.ev2erg(dEg);
        // for (int i = 0; i < a.intervals.size(); i++) {
        // MatterInterval inl = a.getInterval(i);
        // System.out.println("k = " + inl.getWaveNumber(energy) +
        // "\tWidth:" +
        // inl.getWidth() + "\txFromStart = "
        // + a.getSubStackWidth(1, i) + "\tDelta = " + a.getDelta(i,
        // energy));
        // }
        // System.out.println(a);
        // System.out.println("Energy range array length: " +
        // a.getEnergyRange().length + " from " + a.getEnergyRange()[0]
        // + " to "
        // + a.getEnergyRange()[a.getEnergyRange().length - 1] +
        // " ergs");
        // System.out.println("Energy range array length: " +
        // a.getEnergyRange().length + " from " + a.getEnergyRange()[0]
        // + " to "
        // + Utils.erg2ev(a.getEnergyRange()[a.getEnergyRange().length -
        // 1]) +
        // " eV");
        // System.out.println("Eigen energy array length: " +
        // a.getEigenEnergy().length);
        // System.out.println("MatterStack width: " + a.getWidth());
        a.printEigenEnergy(45);
        System.out.println(a.getPropagationMatrix(2, Utils.ev2erg(0.17)));
        System.out.println(a.getWavefunctionCoeffs(2, Utils.ev2erg(0.17))[0] + "    "
                + a.getWavefunctionCoeffs(2, Utils.ev2erg(0.17))[1]);
        System.out.println(a.getPropagationMatrix(0, Utils.ev2erg(0.17)));
        System.out.println(a.getWavefunctionCoeffs(0, Utils.ev2erg(0.17))[0] + "    "
                + a.getWavefunctionCoeffs(0, Utils.ev2erg(0.17))[1]);
        a.plotPropagationCurves();
    }
}
