public class Solution {
    private static double calcF1(double x, double y) {
        return Math.pow(x, 2) / 4 + Math.pow(y, 2) / 9 - 1;
    }

    private static double calcF2(double x, double y) {
        return x - Math.pow(y, 2);
    }

    private static double calcF1DerivX(double x, double y) {
        return x / 2;
    }

    private static double calcF1DerivY(double x, double y) {
        return 2 * y / 9;
    }

    private static double calcF2DerivX(double x, double y) {
        return 1;
    }

    private static double calcF2DerivY(double x, double y) {
        return -2 * y;
    }

    private static double calcFi(double x, double y) {
        return Math.pow(Math.pow(x, 2) / 4 + Math.pow(y, 2) / 9 - 1, 2) + Math.pow(x - Math.pow(y, 2), 2);
    }

    private static double calcFiDerivX(double x, double y) {
        return x * (Math.pow(x, 2) / 4 + Math.pow(y, 2) / 9 - 1) + 2 * (x - Math.pow(y, 2));
    }

    private static double calcFiDerivY(double x, double y) {
        return (4 / 9) * y * (Math.pow(x, 2) / 4 + Math.pow(y, 2) / 9 - 1) - 4 * y * (x - Math.pow(y, 2));
    }

    private static double calcFiDerivT(double x, double y, double t) {
        double fiDerivX = calcFiDerivX(x, y);
        double fiDerivY = calcFiDerivY(x, y);
        return 2 * ((-fiDerivX * (x - t * fiDerivX) / 2 - (2 / 9) * fiDerivY * (y - t * fiDerivY)) *
                (Math.pow(x - t * fiDerivX, 2) / 4 + Math.pow(y - t * fiDerivY, 2) / 9 - 1) +
                (2 * fiDerivY * (y - t * fiDerivY) - fiDerivX) * (-t * fiDerivX - Math.pow(y - t * fiDerivY, 2) + x));
    }

    private static double calcFiDerivTDeriv(double x, double y, double t) {
        double fiDerivX = calcFiDerivX(x, y);
        double fiDerivY = calcFiDerivY(x, y);
        return 2 * ((Math.pow(fiDerivX, 2) / 2 + 2 * Math.pow(fiDerivY, 2) / 9) * (Math.pow(x - t * fiDerivX, 2) / 4 + Math.pow(y - t * fiDerivY, 2) / 9 - 1) -
                4 * Math.pow(fiDerivY, 2) * (x - t * fiDerivX - Math.pow(y - t * fiDerivY, 2)) + 2 * Math.pow(-fiDerivX * (x - t * fiDerivX) / 2 -
                (2 / 9) * fiDerivY * (y - t * fiDerivY), 2) + 2 * Math.pow(2 * fiDerivY * (y - t * fiDerivY) - fiDerivX, 2));
    }

    private static double calcApproxT(double x, double y) {
        double leftBorder = -5;
        double rightBorder = 10;
        double tmp;
        while (rightBorder - leftBorder >= epsilon) {
            tmp = (leftBorder + rightBorder) / 2;
            if (calcFiDerivT(x, y, leftBorder) * calcFiDerivT(x, y, tmp) < 0) {
                rightBorder = tmp;
            } else if (calcFiDerivT(x, y, leftBorder) * calcFiDerivT(x, y, tmp) > 0){
                leftBorder = tmp;
            } else {
                break;
            }
        }
        return (leftBorder + rightBorder) / 2;
    }

    private static double calcT(double x, double y) {
        double tcur;
        double tnext = calcApproxT(x, y);
        do {
            tcur = tnext;
            tnext = tcur - calcFiDerivT(x, y, tcur) / calcFiDerivTDeriv(x, y, tcur);
        } while (Math.abs(tnext - tcur) >= epsilon);
        return tnext;
    }

    private static double[] calcApproxXY() {
        double xcur;
        double xnext = 1;
        double ycur;
        double ynext = 1;
        double tcur;
        double[] result = new double[2];
        do {
            xcur = xnext;
            ycur = ynext;
            tcur = calcT(xcur, ycur);
            xnext = xcur - tcur * calcFiDerivX(xcur, ycur);
            ynext = ycur - tcur * calcFiDerivY(xcur, ycur);
        } while (Math.max(Math.abs(xnext - xcur), Math.abs(ynext - ycur)) >= 0.5);
        result[0] = xnext;
        result[1] = ynext;
        return result;
    }

    private static double calcMatrixNorm(double[][] mtr) {
        double result = -1;
        double sum;
        for(int i = 0; i < mtr.length; i++) {
            sum = 0;
            for (int j = 0; j < mtr[i].length; j++) {
                sum += Math.abs(mtr[i][j]);
            }
            result = Math.max(result, sum);
        }
        return result;
    }

    private static double calcVectorNorm(double[] vector) {
        double result = Math.abs(vector[0]);
        for(int i = 1; i < vector.length; i++) {
            result = Math.max(result, vector[i]);
        }
        return result;
    }

    private static double[] calcVectorSubtract(double[] left, double[] right) {
        double[] result = left.clone();
        for (int i = 0; i < left.length; i++) {
            result[i] -= right[i];
        }
        return result;
    }

    private static double[] calcMatrixMultiplyVector(double[][] mtr, double[] vector) {
        double[] result = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            result[i] = 0;
            for (int j = 0; j < mtr[i].length; j++) {
                result[i] += mtr[i][j] * vector[j];
            }
        }
        return result;
    }

    private static double[][] calcJacobian(double x, double y) {
        double[][] result = new double[2][2];
        result[0][0] = calcF1DerivX(x, y);
        result[0][1] = calcF1DerivY(x, y);
        result[1][0] = calcF2DerivX(x, y);
        result[1][1] = calcF2DerivY(x, y);
        return result;
    }

    private static double[][] calcInverseJacobian(double x, double y) {
        double[][] result = new double[2][2];
        double c = 1 / (calcF1DerivX(x, y) * calcF2DerivY(x, y) - calcF1DerivY(x, y) * calcF2DerivX(x, y));
        result[0][0] = c * calcF2DerivY(x, y);
        result[0][1] = -c * calcF1DerivY(x, y);
        result[1][0] = -c * calcF2DerivX(x, y);
        result[1][1] = c * calcF1DerivX(x, y);
        return result;
    }

    private static void checkTheoremConditions(double x, double y) {
        double[][] jacobian;
        double[][] jacobianInverse;
        double[] fVect = new double[2];
        double[] exactSolution = new double[2];
        double[] approxSolution = new double[2];
        double a1;
        double a2;
        double delta;
        exactSolution[0] = 1.7900855862527592;
        exactSolution[1] = 1.3379408007280291;
        approxSolution[0] = x;
        approxSolution[1] = y;
        delta = calcVectorNorm(calcVectorSubtract(exactSolution, approxSolution));
        fVect[0] = calcF1(exactSolution[0], exactSolution[1]) - calcF1(x, y);
        fVect[1] = calcF2(exactSolution[0], exactSolution[1]) - calcF2(x, y);
        jacobian = calcJacobian(exactSolution[0], exactSolution[1]);
        jacobianInverse = calcInverseJacobian(exactSolution[0], exactSolution[1]);
        a1 = calcMatrixNorm(jacobianInverse);
        a2 = calcVectorNorm(calcVectorSubtract(fVect, calcMatrixMultiplyVector(jacobian, calcVectorSubtract(exactSolution, approxSolution)))) /
                Math.pow(calcVectorNorm(calcVectorSubtract(exactSolution, approxSolution)), 2);
        System.out.println("a1: " + a1);
        System.out.println("a2: " + a2);
        System.out.println("c: " + a1 * a2);
        System.out.println("delta: " + delta);
        System.out.println("b: " + Math.min(1 / (a1 * a2), delta));
    }

    private static double epsilon = 0.00001;

    public static void main(String[] args) {
        double[] xy = calcApproxXY();
        double xcur;
        double xExternal;
        double xnext = xy[0];
        double ycur;
        double yExternal;
        double ynext = xy[1];
        int iterations = 0;
        checkTheoremConditions(xnext, ynext);
        System.out.println("X приближённое: " + xnext);
        System.out.println("Y приближённое: " + ynext);
        do {
            xExternal = xnext;
            yExternal = ynext;
            do {
                xcur = xnext;
                ycur = ynext;
                xnext = xcur - calcF1(xcur, yExternal) / calcF1DerivX(xcur, yExternal);
                ynext = ycur - calcF2(xcur, ycur) / calcF2DerivY(xcur, ycur);
            } while (Math.max(Math.abs(xnext - xcur), Math.abs(ynext - ycur)) >= epsilon / 2);
            iterations++;
        } while (Math.max(Math.abs(xnext - xExternal), Math.abs(ynext - yExternal)) >= epsilon);
        System.out.println("X найденное: " + xnext);
        System.out.println("Y найденное: " + ynext);
        System.out.println("Невязка: " + Math.max(Math.abs(calcF1(xnext, ynext)), Math.abs(calcF2(xnext, ynext))));
        System.out.println("Количество итераций: " + iterations);
    }
}
