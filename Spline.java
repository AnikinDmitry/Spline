//Approximation of exp(x) on [0, 1] by spline

import java.io.FileWriter;
import java.util.Formatter;

public class Spline {
    /*
    sn - spline number:
    spline №1 (ideal spline) has second derivatives at the ends like a test function
    spline №2 (natural spline) has second derivatives at the ends equaled zero
    spline №3 (Lagrange spline) has second derivatives at the ends like Lagrange 
    polynomial of 3 degrees, constructed at 4 edge points
    spline №4 (exotic spline) has continuous third derivatives at the second 
    and penultimate points
    */
    public static void main(String[] args) throws Exception {
        FileWriter out = new FileWriter("out.txt");
        Formatter s = new Formatter();
        int n;
        s.format("Ideal spline:\n");
        for (int i = 0; i < 4; i++) {
            n = (int) (10 * Math.pow(2, i));
            s.format("%f\t%.15f\t%.15f\n", 1.0 / n, splineError(1, n), splineConvergence(1, n));
        }
        s.format("Natural spline:\n");
        for (int i = 0; i < 4; i++) {
            n = (int) (10 * Math.pow(2, i));
            s.format("%f\t%.15f\t%.15f\n", 1.0 / n, splineError(2, n), splineConvergence(2, n));\
        }
        s.format("Lagrange spline:\n");
        for (int i = 0; i < 4; i++) {
            n = (int) (10 * Math.pow(2, i));
            s.format("%f\t%.15f\t%.15f\n", 1.0 / n, splineError(3, n), splineConvergence(3, n));
        }
        s.format("Exotic spline:\n");
        for (int i = 0; i < 4; i++) {
            n = (int) (10 * Math.pow(2, i));
            s.format("%f\t%.15f\t%.15f\n", 1.0 / n, splineError(4, n), splineConvergence(4, n));
        }
        out.write(String.valueOf(s));
        out.close();
        //For work with Gnuplot
        double x;
        FileWriter pointsIdealSpline = new FileWriter("pointsIdealSpline.txt");
        for (int i = 0; i <= 1000; i++) {
            x = (double) i / 1000;
            pointsIdealSpline.write(x + " " + Math.abs(Math.exp(x) - spline(1, 10, x)) + "\n");
        }
        pointsIdealSpline.close();
        FileWriter pointsExoticSpline = new FileWriter("pointsExoticSpline.txt");
        for (int i = 0; i <= 1000; i++) {
            x = (double) i / 1000;
            pointsExoticSpline.write(x + " " + Math.abs(Math.exp(x) - spline(4, 10, x)) + "\n");
        }
        pointsExoticSpline.close();
    }

    public static double spline(int sn, int n, double x) {
        if (sn >= 1 && sn <= 4) {
            double S = 0, h = 1.0 / n;
            double[][] A;
            double[] f;
            double[] M = new double[n + 1];
            if (sn != 4) {
                A = new double[n + 1][n + 1];
                f = new double[n + 1];
                for (int i = 1; i <= n - 1; i++) {
                    A[i - 1][i] = 1;
                    A[i][i] = 4;
                    A[i][i - 1] = 1;
                    f[i] = 6 * (Math.exp((i + 1) * h) + Math.exp((i - 1) * h) - 
                        2 * Math.exp(i * h)) / (h * h);
                }
                A[0][1] = 0;
                A[n - 1][n] = 1;
                A[0][0] = 1;
                A[n][n] = 1;
                if (sn == 1) {
                    f[0] = 1;
                    f[n] = Math.exp(1);
                }
                if (sn == 3) {
                    f[0] = (2 - 5 * Math.exp(h) + 4 * Math.exp(2 * h) - 
                        Math.exp(3 * h)) / (h * h);
                    f[n] = (2 * Math.exp(1) - 5 * Math.exp((n - 1) * h) 
                        + 4 * Math.exp((n - 2) * h) - Math.exp((n - 3) * h)) / (h * h);
                }
                System.arraycopy(sweep(n + 1, A, f), 0, M, 0, n + 1);
            } else {
                A = new double[n - 1][n - 1];
                f = new double[n - 1];
                for (int i = 1; i <= n - 3; i++) {
                    A[i - 1][i] = 1;
                    A[i][i] = 4;
                    A[i][i - 1] = 1;
                    f[i] = 6 * (Math.exp((i + 2) * h) + Math.exp(i * h) - 
                        2 * Math.exp((i + 1) * h)) / (h * h);
                }
                A[0][1] = 0;
                A[n - 3][n - 2] = 1;
                A[0][0] = 1;
                A[n - 2][n - 2] = 1;
                f[0] = (Math.exp(2 * h) + 1 - 2 * Math.exp(h)) / (h * h);
                f[n - 2] = (Math.exp(1) + Math.exp((n - 2) * h) - 
                    2 * Math.exp((n - 1) * h)) / (h * h);
                System.arraycopy(sweep(n - 1, A, f), 0, M, 1, n - 1);
                M[0] = 2 * M[1] - M[2];
                M[n] = 2 * M[n - 1] - M[n - 2];
            }
            for (int i = 1; i <= n; i++)
                if (x >= (i - 1) * h && x <= i * h)
                    S = M[i - 1] * Math.pow(i * h - x, 3) / (6 * h) + M[i] * Math.pow(x - 
                        (i - 1) * h, 3) / (6 * h) + (Math.exp((i - 1) * h) - 
                        M[i - 1] * h * h / 6) * (i * h - x) / h + 
                        (Math.exp(i * h) - M[i] * h * h / 6) * (x - (i - 1) * h) / h;
            return S;
        } else throw new ArithmeticException("Spline isn't defined");
    }

    public static double splineError(int sn, int n) {
        double e = 0;
        double h = 1.0 / (10 * n);
        for (int i = 1; i <= 10 * n; ++i) {
            double le = Math.abs(Math.exp(i * h) - spline(sn, n, i * h));
            if (le > e) e = le;
        }
        return e;
    }

    public static double splineConvergence(int sn, int n) {
        return Math.log(Math.abs((splineError(sn, n) - splineError(sn, 2 * n)) /
                (splineError(sn, 2 * n) - splineError(sn, 4 * n)))) / Math.log(2);
    }

    public static double[] sweep(int N, double[][] A, double[] f) {
        double[] a = new double[N];
        double[] b = new double[N];
        a[0] = -A[0][1] / A[0][0];
        b[0] = f[0] / A[0][0];
        for (int i = 1; i < N - 1; i++) {
            a[i] = -A[i][i + 1] / (A[i][i] + a[i - 1] * A[i][i - 1]);
            b[i] = (f[i] - b[i - 1] * A[i][i - 1]) / (A[i][i] + a[i - 1] * A[i][i - 1]);
        }
        double[] X = new double[N];
        X[N - 1] = (f[N - 1] - b[N - 2] * A[N - 1][N - 2]) / (A[N - 1][N - 1] + 
            a[N - 2] * A[N - 1][N - 2]);
        for (int i = N - 2; i >= 0; i--)
            X[i] = a[i] * X[i + 1] + b[i];
        return X;
    }
}