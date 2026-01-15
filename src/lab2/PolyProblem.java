package lab2;

import org.jgap.*;
import org.jgap.impl.*;

public class PolyProblem {
    private static final int NUM_EVOLUTIONS = 500;
    private static final int NR_COEFICIENTI = 6;
    private static final int GENE_PER_COEF = 16;

    public static void main(String[] args) throws InvalidConfigurationException {
        int totalGenes = NR_COEFICIENTI * GENE_PER_COEF;

        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DeltaFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);

        conf.setFitnessFunction(new PolyFitnessFunction());

        Gene[] sampleGenes = new Gene[totalGenes];
        for (int i = 0; i < totalGenes; i++) {
            sampleGenes[i] = new BooleanGene(conf);
        }

        IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
        conf.setSampleChromosome(sampleChromosome);

        conf.setPopulationSize(1000);

        Genotype population = Genotype.randomInitialGenotype(conf);

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();

            if (i % 50 == 0 || i == NUM_EVOLUTIONS - 1) {
                IChromosome best = population.getFittestChromosome();
                System.out.println("Generatia: " + i + " | Eroare: " + best.getFitnessValue());
            }
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Timp executie: " + (endTime - startTime) + " ms");

        IChromosome bestSolution = population.getFittestChromosome();
        System.out.println("\nCea mai buna solutie gasita:");
        displayIndividual(bestSolution);
    }

    public static void displayIndividual(IChromosome chromosome) {
        double[] coefs = decodeCoefficients(chromosome);

        System.out.println("Polinomul p(x):");
        for (int i = 5; i >= 0; i--) {
            System.out.print((coefs[i] >= 0 ? " + " : " - ") + Math.abs(coefs[i]) + " * x^" + i);
        }
        System.out.println();
        System.out.println("Eroare totala: " + chromosome.getFitnessValue());
    }

    public static double[] decodeCoefficients(IChromosome chromosome) {
        double[] coefs = new double[NR_COEFICIENTI];
        int geneIndex = 0;

        for (int i = 0; i < NR_COEFICIENTI; i++) {
            double base10 = 0;
            for (int k = 0; k < GENE_PER_COEF; k++) {
                boolean bit = (Boolean) chromosome.getGene(geneIndex).getAllele();
                if (bit) {
                    base10 += Math.pow(2, k);
                }
                geneIndex++;
            }
            coefs[i] = (base10 / (Math.pow(2, GENE_PER_COEF) - 1)) * 20 - 10;
        }
        return coefs;
    }

    public static class PolyFitnessFunction extends FitnessFunction {

        public double evaluate(IChromosome chr) {
            double[] coefs = decodeCoefficients(chr);
            double errorSum = 0;
            int steps = 1000;

            for (int k = 0; k <= steps; k++) {
                double x = (k * Math.PI) / steps;
                double polyVal = calculateVal(coefs, x);
                double targetVal = Math.sin(x);

                errorSum += Math.abs(polyVal - targetVal);
            }

            return errorSum;
        }

        private double calculateVal(double[] coefs, double x) {
            double val = 0;
            for (int i = 0; i < coefs.length; i++) {
                val += coefs[i] * Math.pow(x, i);
            }
            return val;
        }
    }
}

//pentru a ajunge la un rezultat bun, a trebuit sa fac niste modifiari
//in cazul initial, doar 50 de indivizi erau prea putini si aveam o eroare de 3700
//am crescut la 1000 de indivizi pentru a putea explora mai multe zone
//am crescut de la 11 la 16 biti pentru a obtine o rezolutie mai buna la coeficientii puterilor mari, astfel am reusit sa trecem sub 500
//cu valorile initiale (50 populatie) e un algoritm rapid dar ramane blocat in erori mari