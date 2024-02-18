package bms;

import beast.base.evolution.tree.IntervalType;

public class BetaSkylineDistributionVersion3 extends AbstractBetaSkylineDistribution {

    public BetaSkylineDistributionVersion3() {
        super();
    }


    // This version updates N after a certain time
    @Override
    public double calculateLogP() {
        logP = 0.0;
        double T = tree.getRoot().getHeight();
        //System.out.print("New line");
        //System.out.print(T);

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);

        double dt_track = 0;                                                           // this variable is first dt after every round and then t_update is deducted until dt_track < t_update
        double t_update = T/populationSizes.getDimension();                            //this variable defines the length of the time intervals
        int counter = 0;                                                               // counter for groups
        double t_rest = 0;                                                             //this variable is needed that we don't use dt_track after the while loop AND again in the next for loop round
        boolean update = false;

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            dt_track += dt;
            //System.out.print(collapsedTreeIntervals.getInterval(i));

            int n = collapsedTreeIntervals.getLineageCount(i);
            t += dt;
            //System.out.print(t + "\t");

            while(dt_track >= t_update && counter < skylinePopulationsInput.get().getDimension()-1){
                logP += -betaCoalescentModel.getTotalCoalRate(n)*(t_update-t_rest)/N;  // here we calculate the logP of the waiting time for the pre-defined interval t_update as long dt_track is bigger than t_update
                dt_track -= t_update;                                                  // in the first round of the while loop we need to update earlier because there is still some dt_track to deduct from the last round
                t_rest = 0;                                                            //now t_rest should be set to 0 again

                counter += 1;
                N = skylinePopulationsInput.get().getValue(counter);
                update = true;
            }

            t_rest = dt_track;

            if (update) {
                logP += -betaCoalescentModel.getTotalCoalRate(n) * t_rest / N;                 // here we calculate the logP for the time interval which is smaller or equal than t_update
                update = false;
            }else{
                logP += -betaCoalescentModel.getTotalCoalRate(n) * dt / N;
            }


            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) - Math.log(N);          // + Binomial.logChoose(n, k)
            }
        }
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    double[] getPopSizes(int gridSize) {
        double[] popSizes = new double[gridSize];
        double T = tree.getRoot().getHeight();

        int group = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int i = 0;
        double t_total = 0;
        double t_update = T/populationSizes.getDimension();
        double dt_track = 0;

        for (int gridIdx = 0; gridIdx < gridSize; gridIdx++) {
            double t = gridIdx * T / (gridSize - 1);

            // TODO: update N as needed
            while((t_total+collapsedTreeIntervals.getInterval(i))<t){

                double dt = collapsedTreeIntervals.getInterval(i);

                // Increment time
                t_total += dt;
                dt_track += dt;
                while(dt_track >= t_update && group < skylinePopulationsInput.get().getDimension()-1){
                    dt_track -= t_update;
                    group += 1;
                    N = skylinePopulationsInput.get().getValue(group);
                }

                if (i == (collapsedTreeIntervals.getIntervalCount()-1)) {  //this if statement ends while loop when we reach the root of the tree
                    break;
                }
                i++;

            }

            popSizes[gridIdx] = N;
        }

        return popSizes;
    }
}