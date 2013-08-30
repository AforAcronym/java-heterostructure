package ru.ioffe.tools;

public interface Function {
	double CALC_TOLERANCE = 1e-16;

	public double call(double c);

	public double call(double a, double b);
}
