package au.edu.rmit.bdm.clustering.trajectory;

import au.edu.rmit.bdm.clustering.mtree.*;
import au.edu.rmit.bdm.clustering.mtree.tests.Data;
import au.edu.rmit.bdm.clustering.mtree.utils.Utils;

public class TrajectoryMtree extends MTree<Data> {
    private static final PromotionFunction<Data> nonRandomPromotion1 = (dataSet, distanceFunction) -> Utils.minMax(dataSet);

    /*
     * initialize our trajectory functions, we can specify the capacity and distance function
     */
    public TrajectoryMtree() {//DEFAULT_MIN_NODE_CAPACITY
        super(20, DistanceFunctions.EBD,
                new ComposedSplitFunction<>(nonRandomPromotion1, new PartitionFunctions.BalancedPartition<>())
        );
    }

    /*
     * initialize our trajectory functions, we can specify the capacity and distance function
     */
    public TrajectoryMtree(int capacity) {//DEFAULT_MIN_NODE_CAPACITY
        super(capacity, DistanceFunctions.EBD,
                new ComposedSplitFunction<>(
                        nonRandomPromotion1,
                        new PartitionFunctions.BalancedPartition<>()
                )
        );
    }

    public void add(Data data) {
        super.add(data);
        _check();
    }

    public boolean remove(Data data) {
        boolean result = super.remove(data);
        _check();
        return result;
    }

    DistanceFunction<? super Data> getDistanceFunction() {
        return distanceFunction;
    }

    public void buildMtree(int[] trajectory, int traid) {
        Data data = new Data(trajectory, traid);
        add(data);
    }

    public void buildHistogram() {
        buildHistogram(this.root);
    }

    public void writeMtree(String folder) {
        writeMtreetoFile(this.root, 1, this.distanceFunction.getID(this.root.getData()), folder);
    }
}
