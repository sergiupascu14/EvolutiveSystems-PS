package lab3;

import java.util.Random;

public class SystemValidation {

    public static void main(String[] args) {
        double[][] A_real = {{0.8, 0.1}, {0.2, 0.7}};
        double[] B_real = {5.0, 3.0};
        double[] C_real = {0.3, -0.7};
        double D_real = 0.0;
        double[] x_real = {2.0, -2.0};

        double[][] A_est = {{0.3586, -0.8058}, {-0.0696, 0.7772}};
        double[] B_est = {-2.0627, 2.2093};
        double[] C_est = {0.4277, -0.1589};
        double D_est = 1.0593;
        double[] x_est = {5.0184, 0.3554};

        int steps = 50;
        System.out.println("--- TEST 1: SEMNAL PERIODIC (SINUS) ---");
        runTest(steps, A_real, B_real, C_real, D_real, x_real.clone(),
                A_est, B_est, C_est, D_est, x_est.clone(), true);

        System.out.println("\n--- TEST 2: SEMNAL ALEATOR (ZGOMOT) ---");
        runTest(steps, A_real, B_real, C_real, D_real, x_real.clone(),
                A_est, B_est, C_est, D_est, x_est.clone(), false);
    }

    public static void runTest(int steps,
                               double[][] Ar, double[] Br, double[] Cr, double Dr, double[] xr,
                               double[][] Ae, double[] Be, double[] Ce, double De, double[] xe,
                               boolean isPeriodic) {

        double totalErrorSq = 0;
        Random rand = new Random();

        System.out.println("Step | Intrarea(u) | Y_Real   | Y_Identificat | Diferenta");
        System.out.println("---------------------------------------------------------");

        for (int k = 0; k < steps; k++) {
            double u;
            if (isPeriodic) {
                u = 2 * Math.sin(k * 0.3);
            } else {
                u = -2 + rand.nextDouble() * 4;
            }

            double y_real = Cr[0] * xr[0] + Cr[1] * xr[1] + Dr * u;
            double xr0_next = Ar[0][0] * xr[0] + Ar[0][1] * xr[1] + Br[0] * u;
            double xr1_next = Ar[1][0] * xr[0] + Ar[1][1] * xr[1] + Br[1] * u;
            xr[0] = xr0_next;
            xr[1] = xr1_next;

            double y_est = Ce[0] * xe[0] + Ce[1] * xe[1] + De * u;
            double xe0_next = Ae[0][0] * xe[0] + Ae[0][1] * xe[1] + Be[0] * u;
            double xe1_next = Ae[1][0] * xe[0] + Ae[1][1] * xe[1] + Be[1] * u;
            xe[0] = xe0_next;
            xe[1] = xe1_next;

            double diff = Math.abs(y_real - y_est);
            totalErrorSq += diff * diff;

            if (k < 15) {
                System.out.printf("%3d  | %8.4f    | %8.4f | %8.4f      | %8.4f%n",
                        k, u, y_real, y_est, diff);
            }
        }

        System.out.println("...");
        System.out.println("Eroare Medie Patratica (MSE) pe " + steps + " pasi: " + (totalErrorSq / steps));
    }
}