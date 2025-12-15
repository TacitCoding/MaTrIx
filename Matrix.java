import java.util.Arrays;

public class Matrix {
    // Class Variables
    public final int numberOfRows;
    public final int numberOfCols;
    public final double[][] values;

    // Provide test cases etc...
    public static void main(String[] args) {
        
        Matrix mat2 = new Matrix(new double[][] {{0, 2, 0}, 
                                                 {3, 0, 0},
                                                 {0, 0, 1}});   
        
        double[][] a = {{1, 2}, {3, 4}, {5, 6}};
        double[][] b = {{6, 7, 6, 7}, {6, 7, 6, 7}};
        //mat2 = byScaler(2, 2, mat2);
        //mat2 = sumOfRows(mat2, 1, 0);
        mat2.printMatrix();
        System.out.println();
        (mat2.inverse()).printMatrix();
    }

    // Constructor
    public Matrix(double[][] matrix) {
        this.values = matrix;
        this.numberOfRows = matrix.length;
        this.numberOfCols = matrix[0].length;
    }

    // Scalar Multiplication
    public Matrix scalarMultiplication(double scaler) {
        double[][] newArr = new double[this.values.length][this.values[0].length];
        for (int row = 0; row < numberOfRows; row++) 
            for (int col = 0; col < numberOfCols; col++)
                newArr[row][col] = this.values[row][col] * scaler;

        return new Matrix(newArr);
    }

    // Addition
    public Matrix add(Matrix other) {
        if (this.numberOfRows != other.numberOfRows || this.numberOfCols != other.numberOfCols) {
            throw new IllegalArgumentException(
                "Matrix dimension mismatch: cannot add a " +
                this.numberOfRows + "x" + this.numberOfRows +
                " matrix with a " +
                other.numberOfCols + "x" + other.numberOfCols +
                " matrix."
            );
        }

        double[][] newArr = new double[this.values.length][this.values[0].length];
        for (int row = 0; row < numberOfRows; row++) 
            for (int col = 0; col < numberOfCols; col++)
                newArr[row][col] = this.values[row][col] + other.values[row][col];

        return new Matrix(newArr);
    }


    public void printMatrix() {
        for (int row = 0; row < this.values.length; row++) {
            for (int col = 0; col < this.values[row].length; col++) {
                System.out.print(this.values[row][col] + " ");
            }
            System.out.println();
        }


    }
    // Subtraction
    public Matrix subtract(Matrix other) {
        if (this.numberOfRows != other.numberOfRows || this.numberOfCols != other.numberOfCols) {
            throw new IllegalArgumentException(
                "Matrix dimension mismatch: cannot subtract a " +
                this.numberOfRows + "by" + this.numberOfRows +
                " matrix with a " +
                other.numberOfCols + "by" + other.numberOfCols +
                " matrix."
            );
        }

        double[][] newArr = new double[this.values.length][this.values[0].length];
        for (int row = 0; row < numberOfRows; row++) 
            for (int col = 0; col < numberOfCols; col++)
                newArr[row][col] = this.values[row][col] - other.values[row][col];

        return new Matrix(newArr);
    }

    // Multiplication
    public Matrix multiply(Matrix other) {
        if (this.values[0].length != other.values.length) {
            throw new IllegalArgumentException(
                "Matrix dimension mismatch: cannot multiply a " +
                this.values.length + "x" + this.values[0].length +
                " matrix with a " +
                other.values.length + "x" + other.values[0].length +
                " matrix."
            );
        }

        double[][] newArr = new double[this.values.length][this.values[0].length];

        for (int row = 0; row < this.numberOfRows; row++) 
            for (int col = 0; col < this.numberOfCols; col++)
                newArr[row][col] = matrixMultiplicationSum(this.values[row], getCollum(other.values, col));

        return new Matrix(newArr);
    }

    public double matrixMultiplicationSum(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++)
            sum += (a[i] * b[i]);

        return sum;
    }

    private double[] getCollum(double[][] a, int collum) {
        double[] newArr = new double[a.length];

        for (int i = 0; i < a.length; i++) 
            newArr[i] = a[i][collum];

        return newArr;
    }

    // Matrices Transposition
    public Matrix transpose() {
        double[][] newArr = new double[this.values[0].length][this.values.length];

        for (int row = 0; row < this.numberOfRows; row++)
            for (int col = 0; col < this.numberOfCols; col++)
                newArr[col][row] = this.values[row][col];

        return new Matrix(newArr);
    }
    
    private static Matrix exchangeRows(Matrix matrix, int row1, int row2) {
        Matrix newMatrix = new Matrix(Arrays.copyOf(matrix.values, matrix.values.length));
        //double[] current2 = matrix.values[row1];
        newMatrix.values[row1] = matrix.values[row2];
        newMatrix.values[row2] = matrix.values[row1];
        return newMatrix;

    }

    private static Matrix byScaler(int row, double scaler, Matrix matrix) {
        Matrix newMatrix = new Matrix(Arrays.copyOf(matrix.values, matrix.values.length));
        for (int i = 0; i < newMatrix.values[row].length; i++) {
            newMatrix.values[row][i] = matrix.values[row][i] * scaler;
        }
        return newMatrix;
        
    }

    private static Matrix sumOfRows(Matrix matrix, int replaceRow, int sumRow, double scaler) {
        Matrix newMatrix = new Matrix(Arrays.copyOf(matrix.values, matrix.values.length));
        for (int i = 0; i < newMatrix.values[replaceRow].length; i++) {
            newMatrix.values[replaceRow][i] += scaler * matrix.values[sumRow][i];
            
        }
        return newMatrix;
        
    }

    private static Matrix produceIdentity(int parem) {
        Matrix identity = new Matrix(new double[parem][parem]);
        for (int row = 0; row < parem; row++) {
            for (int col = 0; col < parem; col++) {
                if (row == col) {identity.values[row][col] = 1.0;}
            }
        }
        return identity;
    }



    public Matrix inverse() {
        if (this.numberOfRows != this.numberOfCols) //checks for equal row and col
            throw new IllegalArgumentException("Non Square Matrix");
        Matrix mat = new Matrix(Arrays.copyOf(this.values, this.values.length)); //mutatable matrix
        Matrix identity = produceIdentity(this.values.length); //identity producing
        for (int col = 0; col < this.numberOfCols; col++) { //this allows me to go to each diaganol based on col, col

            if (mat.values[col][col] == 0.0) { //swap rows if zero
                for (int i = col + 1; i < mat.values.length; i++) {
                    if (mat.values[i][col] != 0) {
                        mat = exchangeRows(mat, col, i);
                        identity = exchangeRows(identity, col, i);
                    }
                }
            }
            if (mat.values[col][col] == 0.0)
                 throw new IllegalArgumentException("Non Invertable");

            identity = byScaler(col, 1.0/mat.values[col][col], identity); //make the col, col 1
            mat = byScaler(col, 1.0/mat.values[col][col], mat);
            for (int rowsBelow = col + 1; rowsBelow < mat.values.length; rowsBelow++) { //go below and do the adding of rows to make everything below zero
                identity = sumOfRows(identity, rowsBelow, col, -1.0 * mat.values[rowsBelow][col]);
                mat = sumOfRows(mat, rowsBelow, col, -1.0 * mat.values[rowsBelow][col]);   
            }
            for (int rowsBelow = col - 1; rowsBelow >= 0; rowsBelow--) { //same thing but now its above
                identity = sumOfRows(identity, rowsBelow, col, -1.0 * mat.values[rowsBelow][col]);
                mat = sumOfRows(mat, rowsBelow, col, -1.0 * mat.values[rowsBelow][col]);
            }
        }
        return identity; 
        
    }

}