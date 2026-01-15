package lab2;

import org.jgap.*;
import org.jgap.impl.*;

public class FuncGrIV {
    private static final int MAX_ALLOWED_EVOLUTIONS = 10;

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

        for (int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++) {
            population.evolve();
            IChromosome bestChrSoFar = population.getFittestChromosome();
            System.out.println("Fitness: " + bestChrSoFar.getFitnessValue() + "," +
                    "individual: " + FitnessFunctionGrIV.Mapping(bestChrSoFar));
        }
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

//se pastreaza cel mai bun individ neschimbat si atunci ajungem repede la convergenta (gen 6-7)
//daca marim populatia convergenta e mai rapida dar creste timpul de calcul pe generatii
//daca marim numarul de generatii, nu va fi mare diferenta dupa o anumita generatie (aprox 8) pentru ca algoritmul atinge un minim local sau global
