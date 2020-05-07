package au.edu.rmit.bdm.clustering.trajectory;

import au.edu.rmit.bdm.clustering.mtree.*;
import au.edu.rmit.bdm.clustering.mtree.tests.Data;
import au.edu.rmit.bdm.clustering.mtree.utils.Utils;

public class TrajectoryMtree extends MTree<Data> {
    private static final PromotionFunction<Data> nonRandomPromotion1 =
            (dataSet, distanceFunction) -> Utils.minMax(dataSet);
    private static final int DEFAULT_MIN_NODE_CAPACITY = 20;

    /*
     * initialize our trajectory functions with DEFAULT_MIN_NODE_CAPACITY
     */
    public TrajectoryMtree() {
        super(DEFAULT_MIN_NODE_CAPACITY, DistanceFunctions.EDGE_BASED_DISTANCE,
                new ComposedSplitFunction<>(nonRandomPromotion1, new PartitionFunctions.BalancedPartition<>()));
    }

    /*
     * initialize our trajectory functions, we can specify the capacity and distance function
     */
    public TrajectoryMtree(int capacity) {
        super(capacity, DistanceFunctions.EDGE_BASED_DISTANCE,
                new ComposedSplitFunction<>(nonRandomPromotion1, new PartitionFunctions.BalancedPartition<>()));
    }

    @Override
    public void add(Data data) {
        super.add(data);
        _check();
    }

    @Override
    public boolean remove(Data data) {
        boolean result = super.remove(data);
        _check();
        return result;
    }

    DistanceFunction<? super Data> getDistanceFunction() {
        return distanceFunction;
    }

    public void buildMtree(int[] trajectory, int trajId) {
        Data data = new Data(trajectory, trajId);
        add(data);
    }

    public void buildHistogram() {
        buildHistogram(this.root);
    }

    public void writeMtree(String folder) {
        writeMtreetoFile(this.root, 1, this.distanceFunction.getID(this.root.getData()), folder);
    }
}
