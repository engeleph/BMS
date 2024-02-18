package bms;

import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import bms.BetaCoalescentModel;
import pitchfork.models.CollapsedTreeIntervals;

public abstract class AbstractBetaSkylineDistribution extends TreeDistribution {

    public Input<CollapsedTreeIntervals> collapsedTreeIntervalsInput = new Input<>(
            "collapsedTreeIntervals",
            "Collapsed tree intervals object.",
            Input.Validate.REQUIRED);

    public Input<BetaCoalescentModel> betaCoalescentModelInput = new Input<>(
            "model",
            "Beta-coalescent model.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> skylinePopulationsInput = new Input<>(
            "skylinePopulations",
            "Parameter containing skyline population sizes.",
            Input.Validate.REQUIRED);

    public Input<IntegerParameter> groupSizesInput = new Input<>(
            "groupSizes",
            "Parameter containing skyline group sizes.",
            Input.Validate.REQUIRED);

    public AbstractBetaSkylineDistribution() {
        treeIntervalsInput.setRule(Input.Validate.FORBIDDEN);
        treeInput.setRule(Input.Validate.OPTIONAL);
    }


    protected CollapsedTreeIntervals collapsedTreeIntervals;
    protected BetaCoalescentModel betaCoalescentModel;
    protected RealParameter populationSizes;
    protected IntegerParameter groupSizes;
    protected Tree tree;

    @Override
    public void initAndValidate() {
        betaCoalescentModel = betaCoalescentModelInput.get();
        collapsedTreeIntervals = collapsedTreeIntervalsInput.get();
        populationSizes = skylinePopulationsInput.get();
        groupSizes = groupSizesInput.get();

        tree = collapsedTreeIntervals.treeInput.get();
        treeInput.setValue(tree, this);

        int groupSizesSum = 0;
        for (int groupSize : groupSizes.getValues())
            groupSizesSum += groupSize;

        if (groupSizesSum != tree.getLeafNodeCount() - 1)
            throw new IllegalArgumentException("Sum of elements of groupSizes input should equal number of leaves - 1.");
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    /**
     * Return array of population sizes at evenly spaced time points.
     * Used for logging.
     *
     * @param gridSize Size of grid on which to compute population size
     * @return array of population sizes
     */
    abstract double[] getPopSizes(int gridSize);

}
