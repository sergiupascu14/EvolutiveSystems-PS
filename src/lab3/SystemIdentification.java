package lab3;

import org.jgap.*;
import org.jgap.impl.*;

public class SystemIdentification {
    private static final int NUM_EVOLUTIONS = 500;
    private static final int NUM_GENES = 11;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DeltaFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);

        conf.setFitnessFunction(new FitnessFunctionIdentification());

        Gene[] sampleGenes = new Gene[NUM_GENES];
        for (int i = 0; i < NUM_GENES; i++) {
            sampleGenes[i] = new DoubleGene(conf, -10, 10);
        }

        IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
        conf.setSampleChromosome(sampleChromosome);

        conf.setPopulationSize(2000);

        Genotype population = Genotype.randomInitialGenotype(conf);

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();

            if (i % 50 == 0) {
                IChromosome best = population.getFittestChromosome();
                System.out.println("Generatia " + i + " | Eroare patratica: " + best.getFitnessValue());
            }
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Timp executie: " + (endTime - startTime) + " ms");

        IChromosome bestSolution = population.getFittestChromosome();
        System.out.println("\n--- REZULTAT IDENTIFICARE ---");
        displayParameters(bestSolution);
    }

    public static void displayParameters(IChromosome chr) {
        Double[] p = new Double[NUM_GENES];
        for(int i=0; i<NUM_GENES; i++) p[i] = (Double) chr.getGene(i).getAllele();

        System.out.println("Sistem Identificat vs Sistem Real (Aproximativ):");
        System.out.println("Matricea A: [" + String.format("%.4f", p[0]) + " " + String.format("%.4f", p[1]) + "]");
        System.out.println("            [" + String.format("%.4f", p[2]) + " " + String.format("%.4f", p[3]) + "]");
        System.out.println("REAL A:     [0.8000 0.1000]");
        System.out.println("            [0.2000 0.7000]");
        System.out.println("------------------------------------------------");
        System.out.println("Matricea B: [" + String.format("%.4f", p[4]) + "]");
        System.out.println("            [" + String.format("%.4f", p[5]) + "]");
        System.out.println("REAL B:     [5.0000]");
        System.out.println("            [3.0000]");
        System.out.println("------------------------------------------------");
        System.out.println("Matricea C: [" + String.format("%.4f", p[6]) + " " + String.format("%.4f", p[7]) + "]");
        System.out.println("REAL C:     [0.3000 -0.7000]");
        System.out.println("------------------------------------------------");
        System.out.println("Scalar D:   " + String.format("%.4f", p[8]));
        System.out.println("REAL D:     0.0000");
        System.out.println("------------------------------------------------");
        System.out.println("Stari Init: x1=" + String.format("%.4f", p[9]) + ", x2=" + String.format("%.4f", p[10]));
        System.out.println("REAL Init:  x1=2.0000, x2=-2.0000");
    }

    public static class FitnessFunctionIdentification extends FitnessFunction {
        private static final int TAU = 20;

        public double evaluate(IChromosome chr) {
            double[] params = new double[NUM_GENES];
            for (int i = 0; i < NUM_GENES; i++) params[i] = (Double) chr.getGene(i).getAllele();

            double errorStep = calculateScenarioError(params, true);

            double errorRamp = calculateScenarioError(params, false);

            return errorStep + errorRamp;
        }

        private double calculateScenarioError(double[] p, boolean isStep) {
            double errorSum = 0;

            double x1_real = 2.0;
            double x2_real = -2.0;

            double x1_est = p[9];
            double x2_est = p[10];

            for (int k = 0; k <= TAU; k++) {
                double u;
                if (isStep) {
                    u = (k <= 5) ? 0 : 1;
                } else {
                    u = 0.25 * k;
                }

                double y_real = 0.3 * x1_real - 0.7 * x2_real;

                double x1_real_next = 0.8 * x1_real + 0.1 * x2_real + 5 * u;
                double x2_real_next = 0.2 * x1_real + 0.7 * x2_real + 3 * u;
                x1_real = x1_real_next;
                x2_real = x2_real_next;

                double y_est = p[6] * x1_est + p[7] * x2_est + p[8] * u;

                double x1_est_next = p[0] * x1_est + p[1] * x2_est + p[4] * u;
                double x2_est_next = p[2] * x1_est + p[3] * x2_est + p[5] * u;

                x1_est = x1_est_next;
                x2_est = x2_est_next;

                errorSum += Math.pow(y_real - y_est, 2);
            }
            return errorSum;
        }
    }
}

//aici e interesant, ori am gresit(dar nu cred), ori am reusit sa gasesc un sistem echivalent, ceea ce e valid pentru ca exista infinitate de sisteme