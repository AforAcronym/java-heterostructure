package ru.ioffe.tools;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public strictfp class ComplexMatrix {

    private Complex[][] matrix;


    public ComplexMatrix(Complex[][] m) throws Exception {
        /*
         * Here (i, j) notation is used. i - a row index j - a column
		 * index j=1 j=2 
		 * i=1 [ [ 11, 12], 
		 * i=2   [ 21, 22] ] 
		 * 
		 *         c o l s 
		 * r  [ [ 11, 12, 13], 
		 * o    [ 21, 22, 23], 
		 * w    [ 31, 32, 33] ] 
		 * s
		 */
        for (int i = 1; i < m.length; i++) {
            if (m[i].length != m[0].length) {
                throw new Exception("The matrix is not consistent!");
            }
        }

        matrix = m;

    }


    /**
     * Create a matrix of nulls
     *
     * @param rowSize size of rows
     * @param colSize size of columns
     * @return
     */
    private static ComplexMatrix emptyMatrix(int rowSize, int colSize) {
        // [ [row], [row], [row], [row]...]
        Complex[][] arr = new Complex[colSize][rowSize];
        try {
            return new ComplexMatrix(arr);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }


    /**
     * Get number of the matrix rows
     *
     * @return
     */
    public int getRowsNum() {
        return matrix.length;
    }


    /**
     * Get number of the matrix columns
     *
     * @return
     */
    public int getColsNum() {
        return matrix[0].length;
    }


    /**
     * Get size of the matrix columns
     *
     * @return
     */
    public int getColSize() {
        return getRowsNum();
    }


    /**
     * Get size of the matrix rows
     *
     * @return
     */
    public int getRowSize() {
        return getColsNum();
    }


    /**
     * Get size of the matrix in a double array of size 2
     *
     * @return [ numberOfCols, numberOfRows]
     */
    public int[] getSize() {
        return new int[]{getColsNum(), getRowsNum()};
    }


    /**
     * Returns a string of size for printing
     *
     * @return Cols x Rows as "NxM"
     */
    public String getSizeString() {
        return "" + this.getColsNum() + "x" + this.getRowsNum();
    }


    /**
     * Get an element value
     *
     * @param ri row index
     * @param ci column index
     * @return
     */
    public Complex getElem(int ri, int ci) {
        try {
            return matrix[ri][ci];
        } catch (Exception e) {
            System.out.println("Invalid indices in getElem(int ri, int ci)");
            e.printStackTrace();
        }
        return null;
    }


    /**
     * Set a particular element
     *
     * @param rowIndex row index
     * @param colIndex column index
     * @param z        a Complex to assign to
     * @throws Exception if passed indexes are naegative or to large
     */
    private void setElem(int rowIndex, int colIndex, Complex z) throws Exception {
        if (rowIndex >= getColsNum() || colIndex >= getRowsNum() || rowIndex < 0 || colIndex < 0) {
            throw new Exception("Invalid indices in setElem(int rowIndex, int colIndex, Complex z)");
        }
        matrix[rowIndex][colIndex] = z;
    }


    /**
     * Generate zero matrix of passed dimensions
     *
     * @param rowSize
     * @param colSize
     * @return
     */
    public static ComplexMatrix zeroMatrix(int rowSize, int colSize) {
        // A row size is number of columns
        // A column size is number of rows
        ComplexMatrix mx = emptyMatrix(rowSize, colSize);

        for (int colIndex = 0; colIndex < rowSize; colIndex++) {
            for (int rowIndex = 0; rowIndex < colSize; rowIndex++) {
                try {
                    mx.setElem(rowIndex, colIndex, new Complex(0, 0));
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        return mx;
    }


    /**
     * Return square n x n matrix with Complex 1+0i on its diagonal and
     * zeros elsewhere
     *
     * @param n
     * @return
     */
    public static ComplexMatrix identityMatrix(int n) {
        ComplexMatrix mx = zeroMatrix(n, n);
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


    /**
     * Check whether we can operate with the current and passed matrix in
     * ways where the same matrices size is needed
     *
     * @param m ComplexMatrix
     * @throws Exception
     */
    private void checkConsistenceAndThrowException(ComplexMatrix m) throws Exception {
        if (!checkConsistence(this, m)) {
            throw new Exception("Matrix dimensions must agree");
        }
    }


    /**
     * Multiply each element in the matrix by a Complex.
     *
     * @param z
     * @return A new ComplexMatrix instance
     * @throws Exception
     */
    public ComplexMatrix multiplyEach(Complex z) {
        ComplexMatrix mx = zeroMatrix(getRowSize(), getColSize());
        try {
            for (int ci = 0; ci < getColSize(); ci++) {
                for (int ri = 0; ri < getRowSize(); ri++) {
                    mx.setElem(ri, ci, this.getElem(ri, ci).times(z));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }

        return mx;
    }


    /**
     * Multiply each element in the matrix by a double
     *
     * @param c double number
     * @return
     * @throws Exception
     */
    public ComplexMatrix multiplyEach(double c) {
        return multiplyEach(Complex.fromDouble(c));
    }


    /**
     * Multiply the matrix (each element) by -1
     *
     * @return -M
     */
    public ComplexMatrix reverseSign() {
        return multiplyEach(-1.0);
    }


    /**
     * Add another matrix to this one. The initial matrix is immutable, and
     * this method returns a new ComplexMatrix object
     *
     * @param b CompleMatrix of the same size
     * @return sum of this matrix and b
     * @throws Exception
     */
    public ComplexMatrix add(ComplexMatrix b) throws Exception {
        checkConsistenceAndThrowException(b);

        int rowSize = getRowSize();
        int colSize = getColSize();
        ComplexMatrix mx = emptyMatrix(rowSize, colSize);

        for (int ri = 0; ri < rowSize; ri++) {
            for (int ci = 0; ci < colSize; ci++) {
                mx.setElem(ri, ci, getElem(ri, ci).add(b.getElem(ri, ci)));
            }
        }
        return mx;
    }


    /**
     * Subtract another matrix from this one
     *
     * @param b CompleMatrix of the same size
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
     * Check on exact dimensions
     *
     * @param rowSize size of a row
     * @param colSize size of a column
     * @return
     */
    public boolean hasDimensions(int rowSize, int colSize) {
        return (getRowSize() == rowSize && getColSize() == colSize);
    }


    /**
     * Divide matrix info four matrices, breaking on specified row and
     * column
     *
     * @param subRowSize rows number of the left-top submatrix
     * @param subColSize columns number of the left-top submatrix
     * @return
     */
    public ComplexMatrix[] divideMatrixIntoFour(int subRowSize, int subColSize) {

        int bigRowSize = getRowSize();
        int bigColSize = getColSize();

		/*
		 * a b c d
		 */
        ComplexMatrix a = emptyMatrix(subRowSize, subColSize);
        ComplexMatrix b = emptyMatrix(bigRowSize - subRowSize, subColSize);
        ComplexMatrix c = emptyMatrix(subRowSize, bigColSize - subColSize);
        ComplexMatrix d = emptyMatrix(bigRowSize - subRowSize, bigColSize - subColSize);

        for (int ri = 0; ri < bigRowSize; ri++) {
            for (int ci = 0; ci < bigColSize; ci++) {

                if (ri < subRowSize && ci < subColSize) {
                    try {
                        a.setElem(ri, ci, this.getElem(ri, ci));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else if (ri >= subRowSize && ci < subColSize) {
                    try {
                        b.setElem(ri - subRowSize, ci, this.getElem(ri, ci));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else if (ri < subRowSize && ci >= subColSize) {
                    try {
                        c.setElem(ri, ci - subColSize, this.getElem(ri, ci));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else {
                    try {
                        d.setElem(ri - subRowSize, ci - subColSize, this.getElem(ri, ci));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        return new ComplexMatrix[]{a, b, c, d};
    }


    /**
     * Glue matrices in one bigger matrix
     *
     * @param a left-top matrix
     * @param b right-top matrix
     * @param c left-bottom matrix
     * @param d right-bottom matrix
     * @return
     */
    private ComplexMatrix joinMatrixFromFour(ComplexMatrix a, ComplexMatrix b, ComplexMatrix c, ComplexMatrix d) {

        int rowSize_a = a.getRowSize();
        int colSize_a = a.getColSize();

        ComplexMatrix mx = emptyMatrix(rowSize_a + b.getRowSize(), colSize_a + c.getColSize());

        for (int ri = 0; ri < mx.getRowSize(); ri++) {
            for (int ci = 0; ci < mx.getColSize(); ci++) {

                if (ri < rowSize_a && ci < colSize_a) {
                    try {
                        mx.setElem(ri, ci, a.getElem(ri, ci));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else if (ri >= rowSize_a && ci < colSize_a) {
                    try {
                        mx.setElem(ri, ci, b.getElem(ri - rowSize_a, ci));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else if (ri < rowSize_a && ci >= colSize_a) {
                    try {
                        mx.setElem(ri, ci, c.getElem(ri, ci - rowSize_a));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else {
                    try {
                        mx.setElem(ri, ci, d.getElem(ri - rowSize_a, ci - colSize_a));
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
            return new ComplexMatrix(new Complex[][]{
                    {
                            a.getElem(0, 0).times(b.getElem(0, 0))
                                    .add(a.getElem(0, 1).times(b.getElem(1, 0))),
                            a.getElem(0, 0).times(b.getElem(0, 1))
                                    .add(a.getElem(0, 1).times(b.getElem(1, 1)))},
                    {
                            a.getElem(1, 0).times(b.getElem(0, 0))
                                    .add(a.getElem(1, 1).times(b.getElem(1, 0))),
                            a.getElem(1, 0).times(b.getElem(0, 1))
                                    .add(a.getElem(1, 1).times(b.getElem(1, 1)))}});
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }


    /**
     * Matrix multiplication. This method chooses optimal algorithm fron
     * implemented ones depending on matrix size
     *
     * @param mx ComplexMatrix of appropriate size for multiplication
     * @return
     * @throws Exception
     */
    public ComplexMatrix multiply(ComplexMatrix mx) throws Exception {
        checkConsistenceAndThrowException(mx);
        // if (getRowsNum() == getColsNum()
        // && 2 * getRowsNum() == (getRowsNum() ^ (getRowsNum() - 1)) +
        // 1
        // && 2 * getColsNum() == (getColsNum() ^ (getColsNum() - 1)) +
        // 1) {
        // return strassenMultiply(mx);
        // } else {
        return straightForwardMultiply(mx);
        // }
    }


    /**
     * NOT IMPLEMENTED Strassen's algorithm for matrices n x n where n =
     * 2^n; as fast as O(n^2.8)
     *
     * @param mx
     * @return
     */
    @SuppressWarnings("unused")
    @Deprecated
    private ComplexMatrix strassenMultiply(ComplexMatrix mx) {
        int n = getRowsNum() / 2; // 0 1 2 3
        ComplexMatrix[] m1 = divideMatrixIntoFour(n, n); // a b c d
        ComplexMatrix[] m2 = mx.divideMatrixIntoFour(n, n); // e f g h

        if (m1[0].getColsNum() == 1) {
            Complex a = m1[0].getElem(0, 0);
            Complex b = m1[1].getElem(0, 0);
            Complex c = m1[2].getElem(0, 0);
            Complex d = m1[3].getElem(0, 0);
            Complex e = m2[0].getElem(0, 0);
            Complex f = m2[1].getElem(0, 0);
            Complex g = m2[2].getElem(0, 0);
            Complex h = m2[3].getElem(0, 0);

            Complex p1 = a.times(f.subtract(h));
            Complex p2 = a.add(b).times(h);
            Complex p3 = c.add(d).times(e);
            Complex p4 = d.times(g.subtract(e));
            Complex p5 = a.add(d).times(e.add(h));
            Complex p6 = b.subtract(d).times(g.add(h));
            Complex p7 = a.subtract(c).times(e.add(f));

            // Complex p1 = m1[0].getElem(0,
            // 0).times(m2[1].getElem(0,
            // 0).subtract(m2[3].getElem(0, 0)));
            // Complex p2 = m1[0].getElem(0, 0).add(m1[1]
            // .getElem(0, 0)).times(m2[3].getElem(0, 0));
            // Complex p3 = m1[2].getElem(0, 0).add(m1[3]
            // .getElem(0, 0)).times(m2[0].getElem(0, 0));
            // Complex p4 = m1[3].getElem(0,
            // 0).times(m2[2].getElem(0,
            // 0).subtract(m2[0].getElem(0, 0)));
            // Complex p5 = m1[0].getElem(0, 0).add(m1[3]
            // .getElem(0, 0)).times(m2[0].getElem(0,
            // 0).add(m2[3].getElem(0, 0)));
            // Complex p6 = m1[1].getElem(0,
            // 0).subtract(m1[3].getElem(0,
            // 0)).times(m2[2].getElem(0, 0).add(m2[3].getElem(0,
            // 0)));
            // Complex p7 = m1[0].getElem(0,
            // 0).subtract(m1[2].getElem(0,
            // 0)).times(m2[0].getElem(0, 0).add(m2[1].getElem(0,
            // 0)));

            // System.out.println(" = = = = = = = = = = = = = = = = ");
            // System.out.println("p1 = " + p1);
            // System.out.println("p2 = " + p2);
            // System.out.println("p3 = " + p3);
            // System.out.println("p4 = " + p4);
            // System.out.println("p5 = " + p5);
            // System.out.println("p6 = " + p6);
            // System.out.println("p7 = " + p7);
            try {
                // return new ComplexMatrix(new Complex[][] {
                // { p5.add(p4).add(p6).add(p2), p1.add(p2) },
                // { p3.add(p4),
                // p1.add(p5).subtract(p3).subtract(p7) } });
                return new ComplexMatrix(new Complex[][]{{p5.add(p4).subtract(p6).subtract(p2), p1.add(p2)},
                        {p3.add(p4), p1.add(p5).subtract(p3).subtract(p7)}});

            } catch (Exception e1) {
                e1.printStackTrace();
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

            // a = p5.add(p4).subtract(p2).add(p6);
            // b = p1.add(p2); // b
            // c = p3.add(p4); // c
            // d = p1.add(p5).subtract(p3).subtract(p7);

            a = p5.add(p4).subtract(p2).subtract(p6);
            b = p1.add(p2); // b
            c = p3.add(p4); // c
            d = p1.add(p5).subtract(p3).subtract(p7);

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        return joinMatrixFromFour(a, b, c, d);
    }


    /**
     * Straight forward matrix multiplication
     *
     * @param mx
     * @return
     */
    private ComplexMatrix straightForwardMultiply(ComplexMatrix mx) {
        ComplexMatrix res = zeroMatrix(getColSize(), getColSize());

        for (int ci = 0; ci < getColSize(); ci++) {
            for (int ci2 = 0; ci2 < getColSize(); ci2++) {
                for (int ri = 0; ri < getRowSize(); ri++) {
                    try {
                        res.setElem(ci2, ci, res.getElem(ci2, ci).add(getElem(ri, ci).times(mx.getElem(ci2, ri))));
                    } catch (Exception e) {
                        e.printStackTrace();
                        return null;
                    }
                }
            }
        }

        return res;
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


    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("ComplexMatrix [\n");
        for (int i = 0; i < getRowsNum(); i++) {
            sb.append("  [");
            for (int j = 0; j < getColsNum(); j++) {
                try {
                    sb.append(getElem(i, j).toString());
                } catch (Exception e) {
                    System.out.println("Exception in ComplexMatrix s.append(getElem(i, j).toString())");
                    e.printStackTrace();
                }
                if (j < getColsNum() - 1) {
                    sb.append(",  ");
                }
            }
            sb.append("]\n");
        }
        sb.append("]");
        return sb.toString();
    }


    public static void main(String[] args) throws Exception {
        ComplexMatrix a = new ComplexMatrix(new Complex[][]{{new Complex(1, 2), new Complex(2, 4)},
                {new Complex(3, 3), new Complex(4, 5)}});
        ComplexMatrix b = new ComplexMatrix(new Complex[][]{{new Complex(1, -3), new Complex(3, 2)},
                {new Complex(2, -1), new Complex(2, 3)}});
        System.out.println(a.getColsNum());
        System.out.println(a.getRowsNum());
        System.out.println(a.getSizeString());
        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("a + b = " + a.add(b));
        System.out.println("a - b = " + a.subtract(b));
        System.out.println("a x b = " + a.multiply(b));

        // Complex one = new Complex(1, 0);
        // Complex zero = new Complex(0, 0);

        ComplexMatrix mxone = identityMatrix(4);

        ComplexMatrix mx = new ComplexMatrix(new Complex[][]{
                {new Complex(1, 1), new Complex(2, 0), new Complex(3, 0), new Complex(4, 0)},
                {new Complex(1, 2), new Complex(1, 0), new Complex(1, 0), new Complex(1, 0)},
                {new Complex(2, 1), new Complex(1, 0), new Complex(1, 0), new Complex(1, 0)},
                {new Complex(2, 2), new Complex(1, 0), new Complex(1, 0), new Complex(1, 0)}});

        // ComplexMatrix mxone = identityMatrix(2);
        //
        // ComplexMatrix mx = new ComplexMatrix(new Complex[][] {
        // { new Complex(2, 0), new Complex(3, 0) },
        // { new Complex(1, 0), new Complex(1, 0) } });

        System.out.println("E  = " + mxone);
        System.out.println("mx = " + mx);
        System.out.println("mx * E = " + mx.multiply(mxone));
        System.out.println("E * mx = " + mxone.multiply(mx));

        System.out.println("mx - mx = " + mx.subtract(mx));
        System.out.println(mx.multiply(mx.multiplyEach(new Complex(2, 0))));

        System.out.println(identityMatrix(3));
    }
}
