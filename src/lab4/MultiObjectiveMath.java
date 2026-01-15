package lab4;

import org.jgap.*;
import org.jgap.impl.*;

public class MultiObjectiveMath {
    private static final int NUM_EVOLUTIONS = 50;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DefaultFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);

        conf.setFitnessFunction(new MOFunction());

        Gene[] sampleGenes = new Gene[1];
        sampleGenes[0] = new DoubleGene(conf, 0, Math.PI / 2);

        IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
        conf.setSampleChromosome(sampleChromosome);

        conf.setPopulationSize(100);

        Genotype population = Genotype.randomInitialGenotype(conf);

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();
        }

        System.out.println("--- Rezultate Finale (Frontul Pareto Aproximat) ---");
        System.out.println("x\t| f1(sin)\t| f2(cos)\t| Fitness (Combinat)");
        System.out.println("----------------------------------------------------");

        for (Object o : population.getPopulation().getChromosomes()) {
            IChromosome chr = (IChromosome) o;
            double x = (Double) chr.getGene(0).getAllele();
            double f1 = Math.sin(x);
            double f2 = Math.cos(x);
            System.out.printf("%.4f\t| %.4f\t| %.4f\t| %.4f%n", x, f1, f2, chr.getFitnessValue());
        }

        IChromosome best = population.getFittestChromosome();
        System.out.println("\nCea mai buna solutie (pentru 2*f1 + 3*f2):");
        double bestX = (Double) best.getGene(0).getAllele();
        System.out.println("x = " + bestX);
        System.out.println("f1 = " + Math.sin(bestX));
        System.out.println("f2 = " + Math.cos(bestX));
    }

    public static class MOFunction extends FitnessFunction {

        public double evaluate(IChromosome chr) {
            double x = (Double) chr.getGene(0).getAllele();

            double f1 = Math.sin(x);
            double f2 = Math.cos(x);

            return 2 * f1 + 3 * f2;
        }
    }
}

//populația tinde sa se aglomereze intr-o anumita zona
//algoritmul va converge catre un punct specific de pe frontul Pareto