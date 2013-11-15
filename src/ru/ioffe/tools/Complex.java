package ru.ioffe.tools;

/**
 * **********************************************************************
 * Compilation: javac Complex.java Execution: java Complex
 * <p/>
 * Data type for complex numbers.
 * <p/>
 * The data type is "immutable" so once you create and initialize a Complex
 * object, you cannot change it. The "final" keyword when declaring re and im
 * enforces this rule, making it a compile-time error to change the .re or .im
 * fields after they've been initialized.
 * <p/>
 * java Complex a = 5.0 + 6.0i b = -3.0 + 4.0i Re(a) = 5.0 Im(a) = 6.0 b + a =
 * 2.0 + 10.0i a - b = 8.0 + 2.0i a * b = -39.0 + 2.0i b * a = -39.0 + 2.0i a /
 * b = 0.36 - 1.52i (a / b) * b = 5.0 + 6.0i conj(a) = 5.0 - 6.0i |a| =
 * 7.810249675906654 tan(a) = -6.685231390246571E-6 + 1.0000103108981198i
 * <p/>
 * ***********************************************************************
 */

public strictfp class Complex {
    private final double re; // the real part
    private final double im; // the imaginary part

    /**
     * Static method to convert a double number (x) to a complex (x + 0i)
     *
     * @return x + 0i
     */
    public static Complex fromDouble(double n) {
        return new Complex(n, 0);
    }

    /**
     * create a new object with the given real and imaginary parts
     *
     * @param real
     * @param imag
     */
    public Complex(double real, double imag) {
        re = real;
        im = imag;
    }

    /**
     * String representation of the invoking Complex object
     */
    public String toString() {
        if (Double.isNaN(im) || Double.isInfinite(im)) {
            return re + " + " + "i" + im;
        }
        if (im == 0) {
            return re + "";
        }
        if (re == 0) {
            return im + "i";
        }
        if (im < 0) {
            return re + " - " + (-im) + "i";
        }
        return re + " + " + im + "i";
    }

    /**
     * Get abs/modulus/magnitude
     *
     * @return
     */
    public double abs() {
        return Math.hypot(re, im);
    }

    /**
     * Get angle/phase/argument
     *
     * @return phase between -pi and pi
     */
    public double phase() {
        return Math.atan2(im, re);
    }

    /**
     * Check on numerical equivalence
     *
     * @param x
     * @return true or false
     */
    public boolean equals(Complex x) {
        return re == x.re && im == x.im;
    }

    /**
     * Check whether real or imaginary part is NaN
     *
     * @return true or false
     */
    public boolean isNan() {
        return Double.isNaN(re) || Double.isNaN(im);
    }

    /**
     * Add to real part, sum of this and argument
     *
     * @param c
     * @return Complex
     */
    public Complex add(double d) {
        return new Complex(re + d, im);
    }

    /**
     * Add argument to this Complex
     *
     * @param c Complex
     * @return Complex
     */
    public Complex add(Complex c) {
        return new Complex(re + c.re, im + c.im);
    }

    /**
     * Subtract argument of type double from real part
     *
     * @param c
     * @return
     */
    public Complex subtract(double c) {
        return new Complex(re - c, im);
    }

    /**
     * Subtract argument
     *
     * @param c
     * @return
     */
    public Complex subtract(Complex c) {
        return new Complex(re - c.re, im - c.im);
    }

    /**
     * Change sign of (a + ib)
     *
     * @return -a - ib
     */
    public Complex reverseSign() {
        return new Complex(-re, -im);
    }

    /**
     * Scalar multiplication on Complex
     *
     * @param b
     * @return
     */
    public Complex times(Complex b) {
        return new Complex(re * b.re - im * b.im, re * b.im + im * b.re);
    }

    /**
     * Scalar multiplication on double
     *
     * @param x
     * @return
     */
    public Complex times(double x) {
        return new Complex(x * re, x * im);
    }

    /**
     * Get a new Complex object whose value is the conjugate of this
     *
     * @return
     */
    public Complex conjugate() {
        return new Complex(re, -im);
    }

    /**
     * Get a new Complex object whose value is the reciprocal of this
     *
     * @return
     */
    public Complex reciprocal() {
        double scale = re * re + im * im;
        return new Complex(re / scale, -im / scale);
    }

    /**
     * Get real part
     *
     * @return
     */
    public double re() {
        return re;
    }

    /**
     * Get imaginary part
     *
     * @return
     */
    public double im() {
        return im;
    }

    /**
     * Returns a / b
     *
     * @param c
     * @return
     */
    public Complex divide(Complex c) {
        return times(c.reciprocal());
    }

    /**
     * Returns a / b
     *
     * @param c
     * @return
     */
    public Complex divide(double c) {
        return new Complex(re / c, im / c);
    }

    /**
     * Complex exponential exp(x)
     *
     * @return
     */
    public Complex exp() {
        return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
    }

    /**
     * Static version of exponent
     */
    public static Complex exp(Complex c) {
        return c.exp();
    }

    /**
     * Complex sine sin(x)
     *
     * @return
     */
    public Complex sin() {
        return new Complex(Math.sin(re) * Math.cosh(im), Math.cos(re) * Math.sinh(im));
    }

    /**
     * Complex cosine cos(x)
     *
     * @return
     */
    public Complex cos() {
        return new Complex(Math.cos(re) * Math.cosh(im), -Math.sin(re) * Math.sinh(im));
    }

    /**
     * Complex tangent tan(x)
     *
     * @return
     */
    public Complex tan() {
        return sin().divide(cos());
    }

    /**
     * A static version of add()
     *
     * @param a Complex
     * @param b Complex
     * @return a + b
     */
    public static Complex add(Complex a, Complex b) {
        return new Complex(a.re + b.re, a.im + b.im);
    }

    /**
     * Construct complex number from its modulus and angle
     *
     * @param abs
     * @param phase in radians
     * @return new Complex object
     */
    public static Complex fromPolar(double abs, double phase) {
        return new Complex(abs * Math.cos(phase), abs * Math.sin(phase));
    }

    /**
     * Power of n
     *
     * @param n
     * @return this^n
     */
    public Complex pow(double n) {
        return fromPolar(Math.pow(abs(), n), phase() * n);
    }

    /**
     * Static version of pow()
     *
     * @param c Complex
     * @param n double
     * @return c^n
     */
    public static Complex pow(Complex c, double n) {
        return fromPolar(Math.pow(c.abs(), n), c.phase() * n);
    }

    /**
     * Static version of pow() for doubles in order to obtain roots of negative
     * number in simplier way
     *
     * @param c
     * @param n
     * @return c^n
     */
    public static Complex pow(double c, double n) {
        Complex result = pow(new Complex(c, 0), n);
        if (c < 0) {
            return new Complex(0.0, result.im);
        }
        return result;
    }

    /**
     * Complex square root
     *
     * @param c
     * @return
     */
    public static Complex sqrt(Complex c) {
        return fromPolar(Math.sqrt(c.abs()), c.phase() / 2);
    }

    /**
     * Complex square root for double
     *
     * @param c
     * @return
     */
    public static Complex sqrt(double c) {
        if (c < 0) {
            return new Complex(0, Math.sqrt(-c));
        }
        return Complex.fromDouble(Math.sqrt(c));
    }

    /**
     * Complex cube root
     *
     * @param c
     * @return
     */
    public static Complex cbrt(Complex c) {
        return fromPolar(Math.cbrt(c.abs()), c.phase() / 3);
    }

    /**
     * Both roots of sqrt(Complex), mathematically correct
     *
     * @param c
     * @return an array of two Complex roots
     */
    public static Complex[] sqrtSet(Complex c) {
        return new Complex[]{fromPolar(Math.sqrt(c.abs()), c.phase() / 2),
                fromPolar(Math.sqrt(c.abs()), c.phase() / 2 + Math.PI)};
    }

    /**
     * Three roots of cbrt(Complex), mathematically correct
     *
     * @param c
     * @return an array of three Complex roots
     */
    public static Complex[] cbrtSet(Complex c) {
        return new Complex[]{fromPolar(Math.cbrt(c.abs()), c.phase() / 3),
                fromPolar(Math.cbrt(c.abs()), (c.phase() + 2 * Math.PI) / 3),
                fromPolar(Math.cbrt(c.abs()), (c.phase() - 2 * Math.PI) / 3)};
    }

    /**
     * Generate polar form of Complex
     *
     * @return
     */
    public String toPolarString() {
        return "" + abs() + " exp(" + phase() + "i)";
    }

    public static void main(String[] args) {
        Complex a = new Complex(1.0, 3.0);
        Complex b = new Complex(-3.0, 4.0);
        Complex apol = Complex.fromPolar(a.abs(), a.phase());
        Complex i = new Complex(0, 1);

        System.out.println("a            = " + a);
        System.out.println("b            = " + b);
        System.out.println("Re(a)        = " + a.re());
        System.out.println("Im(a)        = " + a.im());
        System.out.println("1/a          = " + a.reciprocal());
        System.out.println("a/a          = " + a.divide(a));
        System.out.println("a/abs(a)     = " + a.divide(a.abs()));
        System.out.println("b + a        = " + b.add(a));
        System.out.println("a - b        = " + a.subtract(b));
        System.out.println("a * b        = " + a.times(b));
        System.out.println("b * a        = " + b.times(a));
        System.out.println("a / b        = " + a.divide(b));
        System.out.println("(a / b) * b  = " + a.divide(b).times(b));
        System.out.println("conj(a)      = " + a.conjugate());
        System.out.println("|a|          = " + a.abs());
        System.out.println("tan(a)       = " + a.tan());
        System.out.println("a from polar = " + apol);
        System.out.println("a*a          = " + a.times(a));
        System.out.println("a^2          = " + Complex.pow(a, 2));
        System.out.println("a to polar   = " + a.toPolarString());
        System.out.println("Nan + 1i     = " + new Complex(Double.NaN, 1));
        System.out.println("1 + Nan*i    = " + new Complex(1, Double.NaN));
        System.out.println("Nan + Nan*i  = " + new Complex(Double.NaN, Double.NaN));
        System.out.println("Inf + Nan*i  = " + new Complex(Double.POSITIVE_INFINITY, Double.NaN));
        System.out.println("Inf + Inf*i  = " + new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY));
        System.out.println("sqrt(-1)     = " + Complex.pow(-1, 0.5) + "   " + Complex.sqrt(-1));
        System.out.println("sqrt(1)      = " + Complex.pow(1, 0.5));
        System.out.println("sqrt(3)      = " + Complex.pow(3, 0.5) + "    " + Complex.sqrt(3));
        System.out.println("sqrt(-3)     = " + Complex.pow(-3, 0.5) + "    " + Complex.sqrt(-3));
        System.out.println("sqrt(i)      = " + Complex.pow(i, 0.5) + "    " + Complex.sqrt(i));
        System.out.println("1.0/(1.0-1.0)= " + 1.0 / (1.0 - 1.0));
        System.out.println("1 + a        = " + a.add(1));
        System.out.println("1 - a        = " + a.reverseSign().add(1));
        System.out.println("-(a - 1)     = " + a.add(-1).reverseSign());


    }

}
