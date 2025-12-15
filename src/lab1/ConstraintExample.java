package lab1;

import org.jgap.*;
import org.jgap.impl.DefaultConfiguration;
import org.jgap.impl.IntegerGene;

public class ConstraintExample {
    private static final int NUM_EVOLUTIONS = 100;

    static void main(String[] args) throws InvalidConfigurationException{
        double targetAmount = 1.84;
        Configuration config = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);

        config.setFitnessEvaluator(new DeltaFitnessEvaluator());
        config.setPreservFittestIndividual(true);
        config.setKeepPopulationSizeConstant(true);

        FitnessFunction fitnessFunction = new SampleFitnessFunction(targetAmount);
        config.setFitnessFunction(fitnessFunction);

        Gene[] sampleGenes = new Gene[4];
        sampleGenes[0] = new IntegerGene(config, 0 , 20);
        sampleGenes[1] = new IntegerGene(config, 1 , 30);
        sampleGenes[2] = new IntegerGene(config, 2 , 50);
        sampleGenes[3] = new IntegerGene(config, 3 , 80);
        IChromosome sampleChromosome = new Chromosome(config, sampleGenes);
        config.setSampleChromosome(sampleChromosome);
        config.setPopulationSize(80);
        Genotype population = Genotype.randomInitialGenotype(config);

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();
            IChromosome bestSolutionSoFar = population.getFittestChromosome();
            DisplayIndividual(bestSolutionSoFar);
        }
    }

    public static void DisplayIndividual(IChromosome chromosome){
        System.out.println("Fitness value: " + chromosome.getFitnessValue());
        System.out.println(", Coins: ");
        for (int i = 0; i < 4; i++) {
            System.out.print(SampleFitnessFunction.getNrCoinsAtGene(chromosome, i) + " ");
        }
        System.out.println(", total change: " + SampleFitnessFunction.amountOfChange(chromosome));
    }
}
