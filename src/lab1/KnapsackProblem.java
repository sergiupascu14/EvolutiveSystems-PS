package lab1;

import org.jgap.*;
import org.jgap.impl.BooleanGene;
import org.jgap.impl.DefaultConfiguration;

public class KnapsackProblem {
    private static final int NUM_EVOLUTIONS = 100;

    static void main(String[] args) throws InvalidConfigurationException {
        int targetVolume = 30;
        int objectCount = 10;

        Configuration config = new DefaultConfiguration();
        config.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        config.setFitnessEvaluator(new DefaultFitnessEvaluator());
        config.setPreservFittestIndividual(true);
        config.setKeepPopulationSizeConstant(false);

        FitnessFunction fitnessFunction = new KnapsackFitnessFunction(targetVolume, objectCount);
        config.setFitnessFunction(fitnessFunction);

        Gene[] sampleGenes = new Gene[objectCount];
        for (int i = 0; i < objectCount; i++) {
            sampleGenes[i] = new BooleanGene(config);
        }

        IChromosome sampleChromosome = new Chromosome(config, sampleGenes);
        config.setSampleChromosome(sampleChromosome);
        config.setPopulationSize(50);

        Genotype population = Genotype.randomInitialGenotype(config);

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();
            IChromosome best = population.getFittestChromosome();
            System.out.println("Generatia: " + i);
            displayIndividual(best, targetVolume);
            System.out.println("------------------------");
        }

        IChromosome bestSolution = population.getFittestChromosome();
        System.out.println("\nCea mai buna solutie gasita dupa " + NUM_EVOLUTIONS + " generatii:");
        displayIndividual(bestSolution, targetVolume);
    }

    public static void displayIndividual(IChromosome chromosome, int targetVolume) {
        double totalVol = 0;
        double totalVal = 0;

        System.out.print("Obiecte selectate: ");
        for (int i = 0; i < chromosome.size(); i++) {
            boolean taken = (Boolean) chromosome.getGene(i).getAllele();
            if (taken) {
                int vol = i + 1;
                int val = (i + 1) * (i + 1);
                totalVol += vol;
                totalVal += val;
                System.out.print((i + 1) + " ");
            }
        }

        System.out.println();
        System.out.println("Fitness: " + chromosome.getFitnessValue());
        System.out.println("Volum total: " + totalVol + "/" + targetVolume);
        System.out.println("Valoare totala: " + totalVal);
    }
}
