package lab2;

import org.jgap.*;
import org.jgap.impl.*;

public class MinimLocalTrap {
    private static final int MAX_ALLOWED_EVOLUTIONS = 100;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DeltaFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setPopulationSize(10);
        conf.setFitnessFunction(new FitnessFunctionGrIV());

        int nrGenes = 12;
        IChromosome sampleChromosome = new Chromosome(conf, new BooleanGene(conf), nrGenes);
        conf.setSampleChromosome(sampleChromosome);

        Population initialPop = new Population(conf);
        while (initialPop.size() < 10) {
            IChromosome chr = new Chromosome(conf, new BooleanGene(conf), nrGenes);
            for(int i=0; i<nrGenes; i++) {
                chr.getGene(i).setAllele(Math.random() > 0.5);
            }

            double val = FitnessFunctionGrIV.Mapping(chr);
            if (val < 0 && val > -1.5) {
                initialPop.addChromosome(chr);
            }
        }

        Genotype population = new Genotype(conf, initialPop);

        System.out.println("Start: Toti indivizii sunt fortati in minimul local (x < 0)");

        boolean escaped = false;

        for (int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++) {
            population.evolve();
            IChromosome best = population.getFittestChromosome();
            double x = FitnessFunctionGrIV.Mapping(best);
            double fit = best.getFitnessValue();

            if (!escaped && x > 0.5) {
                System.out.println(">>> A EVADAT din minimul local la generatia: " + i + " <<<");
                escaped = true;
            }

            System.out.println("Gen " + i + " | x = " + String.format("%.4f", x) + " | Fitness: " + String.format("%.4f", fit));
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

//codul are rata de mutatie mare by default (1/12) si atunci evadeaza rapid din minimul local
//am introdus artifical o capcana sa mentinem in minimul local, astfel la start toti indivizii sunt la minimul local
//avem populatie mica pentru a limita algoritmul in explorarea zonelore noi

//cand ramane blocat mai multe generatii, defapt combina doi parinti negativi si atunci si copilul va fi negativ
//reuseste sa evadeze din cauza operatorului de mutatie, da flip la un bit si un numar devine poziti, apoi prin combinare iesim din zona de minim
//apoi, pentru ca acest individ o sa aiba un fitness value mult mai bun si atunci o sa fie pastrat acest individ e pastrat automat