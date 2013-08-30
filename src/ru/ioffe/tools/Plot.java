package ru.ioffe.tools;

import javax.swing.JFrame;
import org.opensourcephysics.frames.PlotFrame;

public class Plot {

	public Plot() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// PlotFrame frame = new PlotFrame("position", "amplitude",
		// "First Plot");
		// frame.setSize(400, 400);
		// for (double x = -10, dx = 0.1; x < 10; x += dx) {
		// frame.append(0, x, Math.sin(x));
		// }
		// frame.setVisible(true);
		// frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		plot(20);
	}

	public static void plot(int size) {
		PlotFrame frame = new PlotFrame("Energy", "Propagation", "Propagation");
		frame.setSize(600, 600);
		for (int i = 0; i < size; i++) {
			// frame.append(0, x[i], y[i]);
			frame.append(0, i, i);
			frame.append(1, i, i * i);
		}

		frame.setConnected(true);
		frame.setMarkerShape(0, 0);
		frame.setMarkerShape(1, 0);

		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	}
}