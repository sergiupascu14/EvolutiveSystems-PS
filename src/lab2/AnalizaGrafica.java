package lab2;

import org.jgap.*;
import org.jgap.impl.*;

public class AnalizaGrafica {
    private static final int MAX_ALLOWED_EVOLUTIONS = 20;
    private static final double TARGET_SOLUTION = 1.10716;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DeltaFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);
        conf.setFitnessFunction(new FitnessFunctionGrIV());

        int nrGenes = 12;
        IChromosome sampleChromosome = new Chromosome(conf, new BooleanGene(conf), nrGenes);
        conf.setSampleChromosome(sampleChromosome);
        conf.setPopulationSize(20);

        Genotype population = Genotype.randomInitialGenotype(conf);

        long totalTime = 0;

        System.out.println("Gen\t| X gasit\t| Abatere (Eroare)\t| Timp (ms)");
        System.out.println("---------------------------------------------------------");

        for (int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++) {
            long startTime = System.nanoTime();
            population.evolve();
            long endTime = System.nanoTime();

            long duration = (endTime - startTime);
            totalTime += duration;

            IChromosome bestChrSoFar = population.getFittestChromosome();
            double valoareX = FitnessFunctionGrIV.Mapping(bestChrSoFar);

            double abatere = Math.abs(valoareX - TARGET_SOLUTION);

            System.out.println(i + "\t| " + String.format("%.5f", valoareX) + "\t| "
                    + String.format("%.5f", abatere) + "\t\t| "
                    + String.format("%.3f", duration / 1000000.0));
        }

        System.out.println("---------------------------------------------------------");
        double avgTime = (double) totalTime / MAX_ALLOWED_EVOLUTIONS / 1000000.0;
        System.out.println("Timp mediu per generatie: " + String.format("%.3f", avgTime) + " ms");
    }

    public static class FitnessFunctionGrIV extends FitnessFunction {
        private final double POSITIVE_BIAS = 5;

        public double evaluate(IChromosome chr) {
            double individ = Mapping(chr);
            double fitness = Math.pow(individ, 4) - 2 * Math.pow(individ, 2) - individ;
            return fitness + POSITIVE_BIAS;
        }

        public static double Mapping(IChromosome chr) {
            double base10 = 0;
            for (int i = 0; i < chr.size(); i++) {
                boolean allele = ((Boolean) chr.getGene(i).getAllele()).booleanValue();
                if (allele) {
                    base10 += Math.pow(2, i);
                }
            }
            double individ = (base10 / (Math.pow(2, 12) - 1)) * 4 - 2;
            return individ;
        }
    }
}

//afisam un tabel usor de preluat in excel pentru a vedea x-ul, abaterea si timp-ul