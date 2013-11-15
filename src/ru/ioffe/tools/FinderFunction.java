package ru.ioffe.tools;

public class FinderFunction {

    public static double findRootByBisection(double x_start, double x_end, Function func) {
        double x_middle = (x_start + x_end) / 2;
        double y_start = func.call(x_start);
        double y_end = func.call(x_end);
        double y_middle = func.call(x_middle);

        if (y_middle == 0.0) {
            return x_middle;
        } else if (y_middle * y_end > 0) {
            x_end = x_middle;
        } else if (y_middle * y_start > 0) {
            x_start = x_middle;
        } else if (y_middle < Function.CALC_TOLERANCE) {
            return x_middle;
        }

        x_middle = (x_start + x_end) / 2;
        y_start = func.call(x_start);
        y_end = func.call(x_end);
        y_middle = func.call(x_middle);
        return x_middle;
    }


    /**
     * @param args
     */
    public static void main(String[] args) {

    }

}
