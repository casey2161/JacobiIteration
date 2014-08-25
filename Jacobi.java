import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Stack;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;


public class Jacobi {
	public static void main(String[] args) throws FileNotFoundException {
		for(int k = 0; k < 10; k++) {
			Matrix values = genMatrix();
			currI = 0;
			currJ = 1;
			Stack<Double> offsets = new Stack<Double>();
			offsets.push(offCalc(values));
			PrintWriter datawriter = new PrintWriter(new File("data" + k + ".txt"));
			values.print(datawriter, 6, 2);
			values.print(6, 2);
			while(offsets.peek() > Math.pow(10, -9)) {
				values = iterate(values);
				//System.out.println(offCalc(values));
				offsets.push(offCalc(values));
				//values.print(6, 2);
			}
			values.print(datawriter, 6, 2);
			while(!offsets.empty()) {
				datawriter.println(offsets.pop());
			}
			datawriter.close();
		}
	}

	private static Matrix iterate(Matrix values) {
		Matrix ret = new Matrix(values.getArrayCopy());
		int[] ind = findMax(ret);
		double[][] tempU = new double[2][2];
		tempU[0][0] = ret.get(ind[0], ind[0]);
		tempU[0][1] = ret.get(ind[1], ind[0]);
		tempU[1][0] = ret.get(ind[0], ind[1]);
		tempU[1][1] = ret.get(ind[1], ind[1]);
		Matrix u = new Matrix(tempU);
		EigenvalueDecomposition evd = u.eig();
		u = evd.getV();
		double[][] givens = new double[5][5];
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				if(i == j) {
					givens[i][j] = 1;
				} else {
					givens[i][j] = 0;
				}
			}
		}
		double mag1 = Math.sqrt(u.get(1, 0) * u.get(1, 0) + u.get(0, 0) * u.get(0, 0));
		givens[ind[0]][ind[0]] = u.get(0, 0) / mag1;
		givens[ind[1]][ind[0]] = u.get(1,  0) / mag1;
		
		mag1 = Math.sqrt(u.get(0, 1) * u.get(0, 1) + u.get(1, 1) * u.get(1, 1));
		givens[ind[1]][ind[1]] = u.get(1, 1) / mag1;
		givens[ind[0]][ind[1]] = u.get(0,  1) / mag1;
		
		Matrix g = new Matrix(givens);
		ret = ret.times(g.transpose());
		ret = g.times(ret);
		return ret;
		
	}
	private static int currI = 0;
	private static int currJ = 1;
	private static int[] findMax(Matrix values) {
		double max = Double.MIN_VALUE;
		int[] ind = new int[2];
	    for(int i = 0; i < 5; i++) {
		    for(int j = i; j < 5; j++) {
		        if(Math.abs(values.get(i, j)) > max) {
		            max = Math.abs(values.get(i, j));
		            ind[0] = i;
		            ind[1] = j;
		        }
		    }
		}
	    return ind;
	}

	private static double offCalc(Matrix values) {
		double sum = 0;
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				if(i != j) {
					sum+= Math.pow(values.get(i, j), 2);
				}
			}
		}
		return sum;
	}

	private static Matrix genMatrix() {
	    Matrix A = Matrix.random(5, 5);
	    for(int i = 0; i < 5; i++) {
	        for(int j = 0; j < 5; j++) {
	            A.set(i, j, A.get(i, j)*100);
	        }
	    }
        // make A symmetric
        for (int m = 0; m < 5; m++) {
            for (int n = m + 1; n < 5; n++) {
                A.set(n, m, A.get(m, n));
            }
        }
        return A;
	}	
	
}
