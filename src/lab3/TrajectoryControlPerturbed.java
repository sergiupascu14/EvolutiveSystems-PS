package lab3;

import org.jgap.*;
import org.jgap.impl.*;
import java.util.Random;

public class TrajectoryControlPerturbed {
    private static final int NUM_EVOLUTIONS = 200;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DeltaFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);

        conf.setFitnessFunction(new FitnessFunctionPerturbed());

        int nrGenes = 13;
        Gene[] sampleGenes = new Gene[nrGenes];
        for (int i = 0; i < nrGenes; i++) {
            sampleGenes[i] = new DoubleGene(conf, -2, 2);
        }

        IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
        conf.setSampleChromosome(sampleChromosome);

        conf.setPopulationSize(1000);

        Genotype population = Genotype.randomInitialGenotype(conf);

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();
            IChromosome best = population.getFittestChromosome();

            if (i % 20 == 0) {
                System.out.println("Generatia " + i + " | Eroare medie (cu zgomot): " + best.getFitnessValue());
            }
        }

        IChromosome bestSolution = population.getFittestChromosome();
        System.out.println("\nSolutia finala (Robusteta la zgomot):");
        displayIndividual(bestSolution);
    }

    public static void displayIndividual(IChromosome chromosome) {
        double[] u = FitnessFunctionPerturbed.Mapping(chromosome);
        double[] y = FitnessFunctionPerturbed.GetY_WithNoise(u);

        System.out.println("Eroare medie estimata: " + chromosome.getFitnessValue());
        System.out.println("\nTraiectoria (o simulare cu zgomot):");
        System.out.println("Step | Control(u) | Iesire(y) | Ref(yref)");
        System.out.println("-----------------------------------------");

        for (int i = 0; i < u.length; i++) {
            System.out.printf("%d    | %8.4f   | %8.4f  | %8.4f%n",
                    i, u[i], y[i], FitnessFunctionPerturbed.getYRef(i));
        }
    }

    public static class FitnessFunctionPerturbed extends FitnessFunction {
        private static final double[] yref = { 0, 0, 2, 2, 1, 1, 1.5, 2, 2, 1, 1, 1, 0, 0 };
        private static final int optimizationHorizon = 13;
        private static final int SIMULATION_RUNS = 20;

        public double evaluate(IChromosome chr) {
            double[] u = Mapping(chr);
            double totalAverageError = 0;

            for(int run=0; run < SIMULATION_RUNS; run++) {
                double[] y = GetY_WithNoise(u);
                double errorSum = 0;
                for (int i = 0; i < optimizationHorizon; i++) {
                    errorSum += Math.abs(y[i] - yref[i]);
                }
                totalAverageError += errorSum;
            }

            return totalAverageError / SIMULATION_RUNS;
        }

        public static double[] Mapping(IChromosome chr) {
            double[] u = new double[optimizationHorizon];
            for (int i = 0; i < chr.size(); i++) {
                u[i] = (Double) chr.getGene(i).getAllele();
            }
            return u;
        }

        public static double[] GetY_WithNoise(double[] u) {
            double[] y = new double[optimizationHorizon];
            double x1 = 0;
            double x2 = 0;
            Random rand = new Random();

            for (int k = 0; k < optimizationHorizon; k++) {
                y[k] = 0.3 * x1 - 0.7 * x2;

                double noise = -0.2 + (0.4 * rand.nextDouble());

                double next_x1 = 0.8 * x1 + 0.1 * x2 + 5 * u[k] + noise;
                double next_x2 = 0.2 * x1 + 0.7 * x2 + 3 * u[k];

                x1 = next_x1;
                x2 = next_x2;
            }
            return y;
        }

        public static double getYRef(int i) {
            return yref[i];
        }
    }
}