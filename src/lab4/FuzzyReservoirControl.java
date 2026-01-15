package lab4;

import org.jgap.*;
import org.jgap.impl.*;

public class FuzzyReservoirControl {
    private static final int NUM_EVOLUTIONS = 200;
    private static final int NUM_GENES = 50;

    public static void main(String[] args) throws InvalidConfigurationException {
        Configuration conf = new DefaultConfiguration();
        Configuration.resetProperty(Configuration.PROPERTY_FITEVAL_INST);
        conf.setFitnessEvaluator(new DefaultFitnessEvaluator());
        conf.setPreservFittestIndividual(true);
        conf.setKeepPopulationSizeConstant(true);

        conf.setFitnessFunction(new FuzzyFitness());

        Gene[] sampleGenes = new Gene[NUM_GENES];
        for (int i = 0; i < NUM_GENES; i++) {
            sampleGenes[i] = new IntegerGene(conf, 1, 5);
        }

        IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
        conf.setSampleChromosome(sampleChromosome);
        conf.setPopulationSize(200);

        Genotype population = Genotype.randomInitialGenotype(conf);

        for (int i = 0; i < NUM_EVOLUTIONS; i++) {
            population.evolve();
            if (i % 50 == 0) {
                IChromosome best = population.getFittestChromosome();
                System.out.println("Gen " + i + " | Fitness J: " + best.getFitnessValue());
            }
        }

        IChromosome best = population.getFittestChromosome();
        System.out.println("\n--- Solutie Optima Identificata ---");
        System.out.println("Fitness Final: " + best.getFitnessValue());

        System.out.println("Tabel Reguli C1 (linii=u, col=h):");
        displayRuleTable(best, 0);

        System.out.println("\nTabel Reguli C2 (linii=u, col=h):");
        displayRuleTable(best, 25);
    }

    private static void displayRuleTable(IChromosome chr, int startIndex) {
        String[] labels = {"NL", "NM", "ZR", "PM", "PL"};
        for (int u = 0; u < 5; u++) {
            for (int h = 0; h < 5; h++) {
                int geneIdx = startIndex + (u * 5 + h);
                Integer val = (Integer) chr.getGene(geneIdx).getAllele();
                System.out.print(labels[val-1] + " ");
            }
            System.out.println();
        }
    }

    public static class FuzzyFitness extends FitnessFunction {
        private final double A_const = 0.1; // Aria
        private final double K1 = 10.0;     // Generator 1 (mai mare)
        private final double K2 = 5.0;      // Generator 2
        private final double Alpha = 1.0;   // constanta debit
        private final double H_REF = 5.0;   // Nivel referinta
        private final int T = 50;           // Orizont timp

        private final double MU = 5.0;

        public double evaluate(IChromosome chr) {
            double totalEnergy = 0;
            double totalHeightError = 0;

            double h = H_REF;

            for (int k = 0; k < T; k++) {
                double u = 3.0 + Math.sin(k * 0.5);

                int u_idx = getFuzzyIndex(u, 0, 6); // u intre 0 si 6
                int h_err_idx = getFuzzyIndex(h - H_REF, -2, 2); // eroare +/- 2

                int ruleIdx = u_idx * 5 + h_err_idx;

                Integer c1_val = (Integer) chr.getGene(ruleIdx).getAllele();
                Integer c2_val = (Integer) chr.getGene(25 + ruleIdx).getAllele();

                double c1 = getCrispValue(c1_val);
                double c2 = getCrispValue(c2_val);

                double s1 = Alpha * c1;
                double s2 = Alpha * c2;

                double p1 = K1 * h * s1;
                double p2 = K2 * h * s2;
                totalEnergy += (p1 + p2);

                double h_next = h + A_const * (u - s1 - s2);

                totalHeightError += Math.abs(h - H_REF);

                h = h_next;
                if (h < 0) h = 0;
            }

            double J = totalEnergy - (MU * totalHeightError);

            if (J < 0) return 1;
            return J;
        }

        private int getFuzzyIndex(double val, double min, double max) {
            double step = (max - min) / 5.0;
            if (val <= min + step) return 0; // NL
            if (val <= min + 2*step) return 1; // NM
            if (val <= min + 3*step) return 2; // ZR
            if (val <= min + 4*step) return 3; // PM
            return 4; // PL
        }

        private double getCrispValue(int fuzzyVal) {
            return (fuzzyVal - 1) * 0.5;
        }
    }
}

//Am creat un controler inteligent pentru vanele barajului
//Cromozomul are 50 de gene cu valori 1-5 (grade fuzzy)
//Fitness-ul masoara doua lucruri opuse: cata energie producem (vrem mult) vs cat gresim nivelul apei (vrem putin)
//Formula: J = Energie - MU * Eroare

//Frontul Pareto e linia care uneste toate solutiile eficiente
//Stanga graficului: Tinem apa stabila, dar nu producem curent
//Dreapta graficului: Producem mult curent, dar golim lacul (eroare mare)
//Nu exista solutie perfecta, algoritmul gaseste variantele de compromis