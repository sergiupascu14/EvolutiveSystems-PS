package lab1;

import org.jgap.FitnessFunction;
import org.jgap.IChromosome;

public class KnapsackFitnessFunction extends FitnessFunction {
    private static int targetVolume;
    private static int objectCount;
    private static final double alpha = 15;

    public KnapsackFitnessFunction(int targetVolume, int objectCount) {
        KnapsackFitnessFunction.targetVolume = targetVolume;
        KnapsackFitnessFunction.objectCount = objectCount;
    }

    @Override
    protected double evaluate(IChromosome chromosome) {
        double totalValue = 0;
        double totalVolume = 0;

        for (int i = 0; i < objectCount; i++) {
            boolean taken = (Boolean) chromosome.getGene(i).getAllele();
            int vol = i + 1;
            int val = (i + 1) * (i + 1);

            if (taken) {
                totalVolume += vol;
                totalValue += val;
            }
        }

        double excess = Math.max(0, totalVolume - targetVolume);
        double penalty = alpha * excess * excess;

        double rawFitness = totalValue - penalty;

        return Math.max(1e-6, rawFitness);
    }
}
