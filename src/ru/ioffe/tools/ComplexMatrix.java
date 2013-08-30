package ru.ioffe.tools;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public strictfp  class ComplexMatrix {

	private Complex[][] matrix;

	public ComplexMatrix(Complex[][] m) throws Exception {
		/*
		 * [ [ 11, 12], [ 21, 22] ]
		 * 
		 * [ [ 11, 12, 13], [ 21, 22, 23], [ 31, 32, 33] ]
		 */
		for (int i = 1; i < m.length; i++) {
			if (m[i].length != m[0].length) {
				throw new Exception("The matrix is not consistent!");
			}
		}

		matrix = m;

	}

	private static ComplexMatrix initEmptyMatrix(int n, int m) {
		Complex[][] arr = new Complex[n][m];
		try {
			return new ComplexMatrix(arr);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	private void setElem(int i, int j, Complex n) throws Exception {
		if (i >= getColsNum() || j >= getRowsNum()) {
			throw new Exception("Invalid indices in setElem(int i, int j, Complex n)");
		}
		matrix[i][j] = n;
	}

	public Complex getElem(int i, int j) {
		try {
			return matrix[i][j];
		} catch (Exception e) {
			System.out.println("Invalid indices in getElem(int i, int j)");
			e.printStackTrace();
		}
		return null;
	}

	public int getRowsNum() {
		return matrix.length;
	}

	public int getColsNum() {
		return matrix[0].length;
	}

	public int[] getSize() {
		return new int[] { this.getColsNum(), this.getRowsNum() };
	}

	public String getSizeString() {
		return "" + this.getColsNum() + "x" + this.getRowsNum();
	}

	public static ComplexMatrix initZeroMatrix(int n, int m) {
		ComplexMatrix mx = initEmptyMatrix(n, m);
		for (int i = 0; i < mx.getRowsNum(); i++) {
			for (int j = 0; j < mx.getRowsNum(); j++) {
				try {
					mx.setElem(i, j, new Complex(0, 0));
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		return mx;
	}

	/**
	 * Return square n x n matrix with Complex ones on its diagonal and zeros
	 * elsewhere
	 * 
	 * @param n
	 * @return
	 */
	public static ComplexMatrix initIdentityMatrix(int n) {
		ComplexMatrix mx = initZeroMatrix(n, n);
		for (int i = 0; i < n; i++) {
			try {
				mx.setElem(i, i, Complex.fromDouble(1.0));
			} catch (Exception e) {
				e.printStackTrace();
				return null;
			}
		}
		return mx;
	}

	/**
	 * Check whether matrices have equal sizes
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean checkConsistence(ComplexMatrix a, ComplexMatrix b) {
		if (a.getSizeString().equals(b.getSizeString())) {
			return true;
		}
		return false;
	}

	private void checkConsistenceAndThrowException(ComplexMatrix m) throws Exception {
		if (!checkConsistence(this, m)) {
			throw new Exception("Matrix dimensions must agree");
		}
	}

	/**
	 * Add another matrix to this one
	 * 
	 * @param b
	 *            CompleMatrix of the same size
	 * @return sum of this matrix and b
	 * @throws Exception
	 */
	public ComplexMatrix add(ComplexMatrix b) throws Exception {
		checkConsistenceAndThrowException(b);
		ComplexMatrix mx = initEmptyMatrix(getColsNum(), getRowsNum());
		for (int i = 0; i < getRowsNum(); i++) {
			for (int j = 0; j < getColsNum(); j++) {
				mx.setElem(i, j, getElem(i, j).add(b.getElem(i, j)));
			}
		}
		return mx;
	}

	/**
	 * Subtract another matrix from this one
	 * 
	 * @param b
	 *            CompleMatrix of the same size
	 * @return sum of this matrix and b
	 * @throws Exception
	 */
	public ComplexMatrix subtract(ComplexMatrix b) {
		try {
			return this.add(b.reverseSign());
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Multiply each element on -1
	 * 
	 * @return
	 */
	public ComplexMatrix reverseSign() {
		ComplexMatrix mx = initEmptyMatrix(getColsNum(), getRowsNum());
		for (int i = 0; i < getRowsNum(); i++) {
			for (int j = 0; j < getColsNum(); j++) {
				try {
					mx.setElem(i, j, getElem(i, j).reverseSign());
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		return mx;
	}

	/**
	 * Check on exact dimensions
	 * 
	 * @param rows
	 *            number of rows
	 * @param cols
	 *            number of columns
	 * @return
	 */
	public boolean hasDimensions(int rows, int cols) {
		return (getRowsNum() == rows && getColsNum() == cols);
	}

	/**
	 * Divide matrix info four matrices, breaking on specified row and column
	 * 
	 * @param rows
	 *            rows number of the left-top matrix
	 * @param cols
	 *            columns number of the left-top matrix
	 * @return
	 */
	public ComplexMatrix[] divideMatrixIntoFour(int rows, int cols) {
		ComplexMatrix a = initEmptyMatrix(rows, cols);
		ComplexMatrix b = initEmptyMatrix(rows, getColsNum() - cols);
		ComplexMatrix c = initEmptyMatrix(getRowsNum() - rows, cols);
		ComplexMatrix d = initEmptyMatrix(getRowsNum() - rows, getColsNum() - cols);
		for (int i = 0; i < getRowsNum(); i++) {
			for (int j = 0; j < getColsNum(); j++) {
				if (i < rows && j < cols) {
					try {
						a.setElem(i, j, this.getElem(i, j));
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else if (i < rows && j >= cols) {
					try {
						b.setElem(i, j - cols, this.getElem(i, j));
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else if (i >= rows && j < cols) {
					try {
						c.setElem(i - rows, j, this.getElem(i, j));
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else {
					try {
						d.setElem(i - rows, j - cols, this.getElem(i, j));
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
		return new ComplexMatrix[] { a, b, c, d };
	}

	/**
	 * Glue matrices in one bigger matrix
	 * 
	 * @param a
	 *            left-top matrix
	 * @param b
	 *            right-top matrix
	 * @param c
	 *            left-bottom matrix
	 * @param d
	 *            right-bottom matrix
	 * @return
	 */
	private ComplexMatrix joinMatrixFromFour(ComplexMatrix a, ComplexMatrix b, ComplexMatrix c, ComplexMatrix d) {
		int rows_a = a.getRowsNum();
		int cols_a = a.getColsNum();

		ComplexMatrix mx = initEmptyMatrix(rows_a + c.getRowsNum(), cols_a + b.getColsNum());

		for (int i = 0; i < mx.getRowsNum(); i++) {
			for (int j = 0; j < mx.getColsNum(); j++) {
				if (i < rows_a && j < cols_a) {
					try {
						mx.setElem(i, j, a.getElem(i, j));
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else if (i < rows_a && j >= cols_a) {
					try {
						mx.setElem(i, j, b.getElem(i, j - rows_a));
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else if (i >= rows_a && j < cols_a) {
					try {
						mx.setElem(i, j, c.getElem(i - rows_a, j));
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else {
					try {
						mx.setElem(i, j, d.getElem(i - rows_a, j - cols_a));
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}

		}
		return mx;
	}

	@Deprecated
	public static ComplexMatrix multyply2x2(ComplexMatrix a, ComplexMatrix b) {
		try {
			return new ComplexMatrix(new Complex[][] {
					{ a.getElem(0, 0).times(b.getElem(0, 0)).add(a.getElem(0, 1).times(b.getElem(1, 0))),
							a.getElem(0, 0).times(b.getElem(0, 1)).add(a.getElem(0, 1).times(b.getElem(1, 1))) },
					{ a.getElem(1, 0).times(b.getElem(0, 0)).add(a.getElem(1, 1).times(b.getElem(1, 0))),
							a.getElem(1, 0).times(b.getElem(0, 1)).add(a.getElem(1, 1).times(b.getElem(1, 1))) } });
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Matrix multiplication. This method chooses optimal algorithm fron
	 * implemented ones depending on matrix size
	 * 
	 * @param mx
	 *            ComplexMatrix of appropriate size for multiplication
	 * @return
	 * @throws Exception
	 */
	public ComplexMatrix multiply(ComplexMatrix mx) throws Exception {
		checkConsistenceAndThrowException(mx);
		if (getRowsNum() == getColsNum() && 2 * getRowsNum() == (getRowsNum() ^ (getRowsNum() - 1)) + 1
				&& 2 * getColsNum() == (getColsNum() ^ (getColsNum() - 1)) + 1) {
			return strassenMultiply(mx);
		} else {
			return straightForwardMultiply(mx);
		}
	}

	/**
	 * Strassen's algorithm for matrices n x n where n = 2^n; as fast as
	 * O(n^2.8)
	 * 
	 * @param mx
	 * @return
	 */
	private ComplexMatrix strassenMultiply(ComplexMatrix mx) {
		int n = getRowsNum() / 2; // 0 1 2 3
		ComplexMatrix[] m1 = divideMatrixIntoFour(n, n); // a b c d
		ComplexMatrix[] m2 = mx.divideMatrixIntoFour(n, n); // e f g h

		if (m1[0].getColsNum() == 1) {
			// a — m1[0].getElem(0, 0)
			// b — m1[1].getElem(0, 0)
			// c — m1[2].getElem(0, 0)
			// d — m1[3].getElem(0, 0)
			// e — m2[0].getElem(0, 0)
			// f — m2[1].getElem(0, 0)
			// g — m2[2].getElem(0, 0)
			// h — m2[3].getElem(0, 0)

			Complex p1 = m1[0].getElem(0, 0).times(m2[1].getElem(0, 0).subtract(m2[3].getElem(0, 0)));
			Complex p2 = m1[0].getElem(0, 0).add(m1[1].getElem(0, 0)).times(m2[3].getElem(0, 0));
			Complex p3 = m1[2].getElem(0, 0).add(m1[3].getElem(0, 0)).times(m2[0].getElem(0, 0));
			Complex p4 = m1[3].getElem(0, 0).times(m2[2].getElem(0, 0).subtract(m2[0].getElem(0, 0)));
			Complex p5 = m1[0].getElem(0, 0).add(m1[3].getElem(0, 0)).times(m2[0].getElem(0, 0).add(m2[3].getElem(0, 0)));
			Complex p6 = m1[1].getElem(0, 0).subtract(m1[3].getElem(0, 0)).times(m2[2].getElem(0, 0).add(m2[3].getElem(0, 0)));
			Complex p7 = m1[0].getElem(0, 0).subtract(m1[2].getElem(0, 0)).times(m2[0].getElem(0, 0).add(m2[1].getElem(0, 0)));

			// System.out.println(" = = = = = = = = = = = = = = = = ");
			// System.out.println("p1 = " + p1);
			// System.out.println("p2 = " + p2);
			// System.out.println("p3 = " + p3);
			// System.out.println("p4 = " + p4);
			// System.out.println("p5 = " + p5);
			// System.out.println("p6 = " + p6);
			// System.out.println("p7 = " + p7);
			try {
				return new ComplexMatrix(new Complex[][] { { p5.add(p4).add(p6).add(p2), p1.add(p2) },
						{ p3.add(p4), p1.add(p5).subtract(p3).subtract(p7) } });
			} catch (Exception e) {
				e.printStackTrace();
				return null;
			}
		}
		ComplexMatrix p1, p2, p3, p4, p5, p6, p7;
		ComplexMatrix a, b, c, d;

		try {
			p1 = m1[0].multiply(m2[1].subtract(m2[3]));
			p2 = m1[0].add(m1[1]).multiply(m2[3]);
			p3 = m1[2].add(m1[3]).multiply(m2[0]);
			p4 = m1[3].multiply(m2[2].subtract(m2[0]));
			p5 = m1[0].add(m1[3]).multiply(m2[0].add(m2[3]));
			p6 = m1[1].subtract(m1[3]).multiply(m2[2].add(m2[3]));
			p7 = m1[0].subtract(m1[2]).multiply(m2[0].add(m2[1]));

			a = p5.add(p4).subtract(p2).add(p6);
			b = p1.add(p2);
			c = p3.add(p4);
			d = p1.add(p5).subtract(p3).subtract(p7);

		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		return joinMatrixFromFour(a, b, c, d);
	}

	/**
	 * Straight forward matrix multiplication (NOT IMPLEMENTED)
	 * 
	 * @param mx
	 * @return
	 */
	private ComplexMatrix straightForwardMultiply(ComplexMatrix mx) {
		// TODO implement
		throw new NotImplementedException();
		// return null;
	}

	/**
	 * Get transposed matrix
	 * 
	 * @return
	 */
	public ComplexMatrix transpose() {
		// TODO implement
		throw new NotImplementedException();
	}

	/**
	 * Multiply each element in the matrix on a Complex
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public ComplexMatrix multiplyEach(Complex c) throws Exception {
		for (int i = 0; i < getRowsNum(); i++) {
			for (int j = 0; j < getColsNum(); j++) {
				this.setElem(i, j, this.getElem(i, j).times(c));
			}
		}
		return this;
	}

	/**
	 * Multiply each element in the matrix on a double
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public ComplexMatrix multiplyEach(double c) throws Exception {
		return multiplyEach(Complex.fromDouble(c));
	}

	public String toString() {
		StringBuffer s = new StringBuffer();
		s.append("ComplexMatrix [\n");
		for (int i = 0; i < getRowsNum(); i++) {
			s.append("  [");
			for (int j = 0; j < getColsNum(); j++) {
				try {
					s.append(getElem(i, j).toString());
				} catch (Exception e) {
					System.out.println("Exception in ComplexMatrix s.append(getElem(i, j).toString())");
					e.printStackTrace();
				}
				if (j < getColsNum() - 1) {
					s.append(",  ");
				}
			}
			s.append("]\n");
		}
		s.append("]");
		return s.toString();
	}

	public static void main(String[] args) throws Exception {
		// ComplexMatrix a = new ComplexMatrix(new Complex[][] {
		// { new Complex(1, 2), new Complex(2, 4) },
		// { new Complex(3, 3), new Complex(4, 5) } });
		// ComplexMatrix b = new ComplexMatrix(new Complex[][] {
		// { new Complex(1, -3), new Complex(3, 2) },
		// { new Complex(2, -1), new Complex(2, 3) } });
		// System.out.println(a.getColsNum());
		// System.out.println(a.getRowsNum());
		// System.out.println(a.getSizeString());
		// System.out.println("a = " + a);
		// System.out.println("b = " + b);
		// System.out.println("a + b = " + a.add(b));
		// System.out.println("a - b = " + a.substract(b));
		// System.out.println("a x b = " + a.multiply(b));
		//
		// Complex one = new Complex(1, 0);
		// Complex zero = new Complex(0, 0);
		// ComplexMatrix mxone = new ComplexMatrix(new Complex[][] { { one,
		// zero, zero, zero }, { zero, one, zero, zero },
		// { zero, zero, one, zero }, { zero, zero, zero, one } });
		// ComplexMatrix mx = new ComplexMatrix(new Complex[][] {
		// { new Complex(1, 1), new Complex(2, 0), new Complex(3, 0), new
		// Complex(4, 0) },
		// { new Complex(1, 2), new Complex(1, 0), new Complex(1, 0), new
		// Complex(1, 0) },
		// { new Complex(2, 1), new Complex(1, 0), new Complex(1, 0), new
		// Complex(1, 0) },
		// { new Complex(2, 2), new Complex(1, 0), new Complex(1, 0), new
		// Complex(1, 0) } });
		// System.out.println("mx = " + mx);
		// System.out.println("mx * E = " + mx.multiply(mxone));
		// System.out.println(mx.substract(mx));
		// System.out.println(mx.multiply(mx.multiplyEach(new Complex(2, 0))));

		System.out.println(initIdentityMatrix(3));
	}
}
