package lab3;

import org.jgap.*;
import org.jgap.impl.*;

public class TrajectoryControlDynamicStop {
    private static final int MAX_ALLOWED_EVOLUTIONS = 500;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DeltaFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);

        conf.setFitnessFunction(new FitnessFunctionTrajectory());

        int nrGenes = 13;
        Gene[] sampleGenes = new Gene[nrGenes];
        for (int i = 0; i < nrGenes; i++) {
            sampleGenes[i] = new DoubleGene(conf, -2, 2);
        }

        IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
        conf.setSampleChromosome(sampleChromosome);

        conf.setPopulationSize(1000);

        Genotype population = Genotype.randomInitialGenotype(conf);

        double previousBestFitness = -1;
        int generationsCount = 0;

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++) {
            population.evolve();
            IChromosome best = population.getFittestChromosome();
            double currentFitness = best.getFitnessValue();

            if (i > 0 && i % 10 == 0) {
                if (previousBestFitness != -1) {
                    double improvement = (previousBestFitness - currentFitness) / previousBestFitness;

                    System.out.println("Gen " + i + " | Fitness: " + currentFitness +
                            " | Imbunatatire vs Gen " + (i-10) + ": " + String.format("%.2f%%", improvement * 100));

                    if (improvement < 0.05) {
                        System.out.println(">>> STOP: Performanta nu a crescut cu >5% in ultimele 10 generatii.");
                        generationsCount = i;
                        break;
                    }
                } else {
                    System.out.println("Gen " + i + " | Fitness: " + currentFitness);
                }

                previousBestFitness = currentFitness;
            }
            generationsCount = i;
        }

        long endTime = System.currentTimeMillis();
        System.out.println("\nAlgoritmul s-a oprit dupa " + generationsCount + " generatii.");
        System.out.println("Timp executie: " + (endTime - startTime) + " ms");

        IChromosome bestSolution = population.getFittestChromosome();
        System.out.println("\nSolutia finala:");
        displayIndividual(bestSolution);
    }

    public static void displayIndividual(IChromosome chromosome) {
        double[] u = FitnessFunctionTrajectory.Mapping(chromosome);
        double[] y = FitnessFunctionTrajectory.GetY(u);

        System.out.println("Eroare totala: " + chromosome.getFitnessValue());
        System.out.println("\nTraiectoria finala:");
        System.out.println("Step | Control(u) | Iesire(y) | Ref(yref)");
        System.out.println("-----------------------------------------");

        for (int i = 0; i < u.length; i++) {
            System.out.printf("%d    | %8.4f   | %8.4f  | %8.4f%n",
                    i, u[i], y[i], FitnessFunctionTrajectory.getYRef(i));
        }
    }

    public static class FitnessFunctionTrajectory extends FitnessFunction {
        private static final double[] yref = { 0, 0, 2, 2, 1, 1, 1.5, 2, 2, 1, 1, 1, 0, 0 };
        private static final int optimizationHorizon = 13;

        public double evaluate(IChromosome chr) {
            double[] u = Mapping(chr);
            double[] y = GetY(u);
            double errorSum = 0;
            for (int i = 0; i < optimizationHorizon; i++) {
                errorSum += Math.abs(y[i] - yref[i]);
            }
            return errorSum;
        }

        public static double[] Mapping(IChromosome chr) {
            double[] u = new double[optimizationHorizon];
            for (int i = 0; i < chr.size(); i++) {
                u[i] = (Double) chr.getGene(i).getAllele();
            }
            return u;
        }

        public static double[] GetY(double[] u) {
            double[] y = new double[optimizationHorizon];
            double x1 = 0; double x2 = 0;
            for (int k = 0; k < optimizationHorizon; k++) {
                y[k] = 0.3 * x1 - 0.7 * x2;
                double next_x1 = 0.8 * x1 + 0.1 * x2 + 5 * u[k];
                double next_x2 = 0.2 * x1 + 0.7 * x2 + 3 * u[k];
                x1 = next_x1; x2 = next_x2;
            }
            return y;
        }

        public static double getYRef(int i) { return yref[i]; }
    }
}