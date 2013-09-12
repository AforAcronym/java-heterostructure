/*
 * 
 */
package ru.ioffe.tools;

/**
 * An Interval class for numeric methods
 * 
 * @author Evgeny Shevchenko
 */
public strictfp class Interval
{

	protected double	x_start;
	protected double	x_end;
	protected double	y;
	protected String	label	= "";



	public String getLabel()
	{
		return label;
	}



	public void setLabel(String label)
	{
		this.label = label;
	}



	public Interval()
	{
		x_start = 0;
		x_end = 1;
		y = 1;

	}



	/**
	 * Initializing Interval object without x axis data (x defaults from 0
	 * to 1)
	 * 
	 * @param y_start
	 * @param angle
	 */
	public Interval(double _xs, double _xe, double _y)
	{
		x_start = _xs;
		x_end = _xe;
		y = _y;
	}



	public final double getWidth()
	{
		return x_end - x_start;
	}



	public void setY(double y)
	{
		this.y = y;
	}



	public double getXStart()
	{
		return x_start;
	}



	public void setXStart(double x_start)
	{
		double width = getWidth();
		this.x_start = x_start;
		this.x_end = x_start + width;
	}



	public void setXEnd(double x_end)
	{
		this.x_end = x_end;
	}



	public double getXEnd()
	{
		return x_end;
	}



	public double getY()
	{
		return y;
	}



	/**
	 * Return coordinate of the middle
	 * 
	 * @return
	 */
	public double getX()
	{
		return (x_start + x_end) / 2;
	}



	/**
	 * Get x coordinates for the curve building of wave functions and
	 * potential profiles. This method doesn't require equality of intervals
	 * in a stack
	 * 
	 * @param density
	 *                number of dots per 1 cm
	 * @return x coordinates in the interval instance for wave functions and
	 *         potential profiles building
	 */
	public double[] getGrid(long density)
	{
		int numberOfDots = (int) (density * getWidth());
		double gridStep = 1 / (double) density;
		double[] grid = new double[numberOfDots];
		for (int i = 0; i < numberOfDots; i++)
		{
			grid[i] = x_start + ((double) i + 0.5) * gridStep;
		}
		return grid;
	}



	/**
	 * Get local x coordinates of the interval for the curve building of
	 * wave functions and potential profiles. This method doesn't require
	 * equality of intervals in a stack
	 * 
	 * @param density
	 *                number of dots per 1 cm
	 * @return x coordinates in the interval instance for wave functions and
	 *         potential profiles building
	 */
	public double[] getLocalGrid(long density)
	{
		int numberOfDots = (int) (density * getWidth());
		double gridStep = 1 / (double) density;
		double[] grid = new double[numberOfDots];
		for (int i = 0; i < numberOfDots; i++)
		{
			grid[i] = ((double) i + 0.5) * gridStep;
		}
		return grid;
	}



	@Override
	public String toString()
	{
		return String.format("X: [%.3f, %.3f]\tY: %.3f", x_start, x_end, y);
	}



	public static void main(String[] args)
	{
		double x0 = 0;
		double step = 0.5;
		double y = 3;

		Interval i1 = new Interval(x0, x0 + step, y);
		Interval i2 = new Interval(i1.getXEnd(), i1.getXEnd() + step, i1.getY() / 2);
		Interval i3 = new Interval(i2.getXEnd(), i2.getXEnd() + step, i2.getY() * 2);

		System.out.println(i1);
		System.out.println(i2);
		System.out.println(i3);
	}
}
