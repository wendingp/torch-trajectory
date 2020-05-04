package au.edu.rmit.bdm.Torch.queryEngine.similarity;

import au.edu.rmit.bdm.Torch.base.helper.GeoUtil;
import au.edu.rmit.bdm.Torch.base.model.TrajEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Comparator;
import java.util.List;

import static au.edu.rmit.bdm.Torch.base.helper.GeoUtil.min;

/**
 * Distances usage: DistanceFunction<TorVertex, TorVertex> distFunc = (p1, p2)
 * -> GeoUtil.distance(p1, p2); Comparator<TorVertex> comparator = (p1, p2) -> {
 * double dist = GeoUtil.distance(p1, p2); if (dist < 8) return 0; return 1; };
 * SimilarityFunction<TorVertex> SIM_MEASURE = new
 * SimilarityFunction<>(distFunc, comparator);
 *
 * @author forrest0402, wendingp
 */
public class SimilarityFunction<T extends TrajEntry> {

    private static Logger logger = LoggerFactory.getLogger(SimilarityFunction.class);
    public static final SimilarityFunction<TrajEntry> DEFAULT;

    static { // TODO ?
        final int EPSILON = 50;
        DistanceFunction<TrajEntry, TrajEntry> distFunc = GeoUtil::distance;
        Comparator<TrajEntry> comparator = (p1, p2) -> {
            if (GeoUtil.distance(p1, p2) <= EPSILON) {
                return 0;
            }
            return 1;
        };

        DEFAULT = new SimilarityFunction<>(distFunc, comparator);
    }

    private final DistanceFunction distFunc;

    public Comparator<T> comparator;

    public SimilarityFunction(DistanceFunction<T, T> distFunc, Comparator<T> comparator) {
        this.distFunc = distFunc;
        this.comparator = comparator;
    }

    /**
     * ED
     *
     * @param traj1
     * @param traj2
     * @return
     */
    public double EuclideanDistance(List<T> traj1, List<T> traj2) { // ED O(mn)
        if (traj1.size() > traj2.size()) {
            List<T> tmp = traj1;
            traj1 = traj2;
            traj2 = tmp;
        }
        int n = traj1.size(), m = traj2.size();
        double minDistance = Double.MAX_VALUE;
        for (int st = 0; st + n <= m; ++st) {
            minDistance = Math.min(minDistance, EuclideanDistanceSameLength(traj1, traj2.subList(st, st + n)));
        }
        return minDistance;
    }

    private double EuclideanDistanceSameLength(final List<T> traj1, final List<T> traj2) { // ED of 2 same size
        // trajectories O(n)
        if (traj1.size() != traj2.size()) {
            throw new IllegalArgumentException("T1 should be of the same length as T2");
        }
        double dist = 0.0;
        for (int i = 0; i < traj1.size(); ++i) {
            dist += Math.sqrt(distFunc.apply(traj1.get(i), traj2.get(i)));
        }
        // return dist; // original
        return dist / traj1.size(); // TODO should it be average here?
    }

    public double LongestCommonSubsequence(List<T> traj1, List<T> traj2, int theta) { // LCSS
        if (traj1 == null || traj2 == null || traj1.size() == 0 || traj2.size() == 0) {
            return 0;
        }

        int[][] dp = new int[traj1.size()][traj2.size()];

        if (comparator.compare(traj1.get(0), traj2.get(0)) == 0) {
            dp[0][0] = 1;
        }

        for (int i = 1; i < traj1.size(); ++i) {
            if (comparator.compare(traj1.get(i), traj2.get(0)) == 0) {
                dp[i][0] = 1;
            } else {
                dp[i][0] = dp[i - 1][0];
            }
        }

        for (int i = 1; i < traj2.size(); ++i) {
            if (comparator.compare(traj2.get(i), traj1.get(0)) == 0) {
                dp[0][i] = 1;
            } else {
                dp[0][i] = dp[0][i - 1];
            }
        }

        for (int i = 1; i < traj1.size(); ++i) {
            for (int j = 1; j < traj2.size(); ++j) {
                if (Math.abs(i - j) <= theta) { // TODO ?
                    if (comparator.compare(traj1.get(i), traj2.get(j)) == 0) {
                        dp[i][j] = 1 + dp[i - 1][j - 1];
                    } else {
                        dp[i][j] = Math.max(dp[i - 1][j], dp[i][j - 1]);
                    }
                }
            }
        }
        return dp[traj1.size() - 1][traj2.size() - 1];
    }

    /**
     * ERP
     *
     * @param traj1 trajectory 1
     * @param traj2 trajectory 2
     * @param g     TODO: ?
     * @return distance
     */
    public double EditDistanceWithRealPenalty(final List<T> traj1, final List<T> traj2, final T g) { // TODO what is g?
        if (traj1 == null || traj1.size() == 0) { // TODO: why allow null?
            double res = 0.0;
            if (traj2 != null) {
                for (T t : traj2) {
                    res += distFunc.apply(t, g);
                }
            }
            return res;
        }
        if (traj2 == null || traj2.size() == 0) {
            double res = 0.0;
            for (T t : traj1) {
                res += distFunc.apply(t, g);
            }
            return res;
        }

        double[][] dp = new double[traj1.size() + 1][traj2.size() + 1];
        for (int i = 1; i <= traj1.size(); ++i) {
            dp[i][0] = distFunc.apply(traj1.get(i - 1), g) + dp[i - 1][0];
        }
        for (int j = 1; j <= traj2.size(); ++j) { //
            dp[0][j] = distFunc.apply(traj2.get(j - 1), g) + dp[0][j - 1];
        }

        for (int i = 1; i <= traj1.size(); ++i) {
            for (int j = 1; j <= traj2.size(); ++j) {
                dp[i][j] = min(dp[i - 1][j - 1] + distFunc.apply(traj1.get(i - 1), traj2.get(j - 1)),
                        dp[i - 1][j] + distFunc.apply(traj1.get(i - 1), g),
                        dp[i][j - 1] + distFunc.apply(g, traj2.get(j - 1)));
            }
        }
        return dp[traj1.size()][traj2.size()] / Math.max(traj1.size(), traj2.size());
    }

    /**
     * EDR
     *
     * @param traj1
     * @param traj2
     * @return
     */
    public double EditDistanceOnRealSequence(final List<T> traj1, final List<T> traj2) {
        if (traj1.isEmpty() || traj2.isEmpty()) {
            if (traj1.isEmpty() && traj2.isEmpty())
                return 0;
            return Math.max(traj1.size(), traj2.size());
        }
        int[][] dp = new int[traj1.size() + 1][traj2.size() + 1];

        for (int i = 1; i <= traj1.size(); ++i) {
            dp[i][0] = i;
        }
        for (int j = 1; j <= traj2.size(); ++j) {
            dp[0][j] = j;
        }

        for (int i = 1; i <= traj1.size(); ++i) {
            for (int j = 1; j <= traj2.size(); ++j) {
                int subCost = 1;
                if (comparator.compare(traj1.get(i - 1), traj2.get(j - 1)) == 0) {
                    subCost = 0;
                }
                dp[i][j] = min(dp[i - 1][j - 1] + subCost, dp[i - 1][j] + 1, dp[i][j - 1] + 1);
            }
        }
        return dp[traj1.size()][traj2.size()];
    }

    /**
     * DTW
     *
     * @param traj1
     * @param traj2
     * @return
     */
    public double DynamicTimeWarping(List<T> traj1, List<T> traj2) {
        if (traj1.size() == 0 && traj2.size() == 0) {
            return 0;
        }
        if (traj1.size() == 0 || traj2.size() == 0) {
            return Integer.MAX_VALUE;
        }

        double[][] dp = new double[traj1.size() + 1][traj2.size() + 1];

        for (int i = 1; i <= traj1.size(); ++i) {
            dp[i][0] = Integer.MAX_VALUE;
        }
        for (int j = 1; j <= traj2.size(); ++j) {
            dp[0][j] = Integer.MAX_VALUE;
        }

        for (int i = 1; i <= traj1.size(); ++i) {
            for (int j = 1; j <= traj2.size(); ++j) {
                dp[i][j] = distFunc.apply(traj1.get(i - 1), traj2.get(j - 1))
                        + min(dp[i - 1][j - 1], dp[i - 1][j], dp[i][j - 1]);
            }
        }
        return dp[traj1.size()][traj2.size()];
    }

    /**
     * Hausdorff
     *
     * @param traj1
     * @param traj2
     * @return
     */
    public double Hausdorff(List<T> traj1, List<T> traj2) { // refactored to O(mn)
        double[][] distMatrix = new double[traj2.size()][traj1.size()]; // TODO why traj2 then traj1 ?
        for (int i = 0; i < distMatrix.length; ++i) {
            for (int j = 0; j < distMatrix[0].length; ++j) {
                distMatrix[i][j] = distFunc.apply(traj1.get(j), traj2.get(i));
            }
        }

        double maxMinDistance = Double.MIN_VALUE;
        for (int i = 0; i < distMatrix.length; ++i) {
            double rowMin = Double.MAX_VALUE;
            for (int j = 0; j < distMatrix[0].length; ++j) {
                rowMin = Math.min(rowMin, distMatrix[i][j]);
            }
            maxMinDistance = Math.max(maxMinDistance, rowMin);
        }
        for (int i = 0; i < distMatrix[0].length; ++i) {
            double colMin = Double.MAX_VALUE;
            for (int j = 0; j < distMatrix.length; ++j) {
                colMin = Math.min(colMin, distMatrix[j][i]);
            }
            maxMinDistance = Math.max(maxMinDistance, colMin);
        }
        return maxMinDistance;
    }

    /**
     * Frechet distance: min length of leash required to connect two separate paths
     * Time complexity: O(mn)
     *
     * @param traj1
     * @param traj2
     * @return
     */
    public double Frechet(final List<T> traj1, final List<T> traj2) {
        double[][] ca = new double[traj2.size()][traj1.size()]; // TODO ?
        for (int i = 0; i < traj2.size(); ++i) {
            for (int j = 0; j < traj1.size(); ++j) {
                ca[i][j] = -1.0D;
            }
        }
        return FrechetHelper(traj2.size() - 1, traj1.size() - 1, ca, traj1, traj2);
    }

    private double FrechetHelper(final int i, final int j, final double[][] ca,
                                 final List<T> traj1, final List<T> traj2) {
        if (ca[i][j] > -1.0D) {
            return ca[i][j];
        }
        if (i == 0 && j == 0) {
            ca[i][j] = distFunc.apply(traj1.get(0), traj2.get(0));
        } else if (j == 0) {
            ca[i][j] = Math.max(FrechetHelper(i - 1, 0, ca, traj1, traj2), distFunc.apply(traj2.get(i), traj1.get(0)));
        } else if (i == 0) {
            ca[i][j] = Math.max(FrechetHelper(0, j - 1, ca, traj1, traj2), distFunc.apply(traj2.get(0), traj1.get(j)));
        } else {
            ca[i][j] = Math.max(distFunc.apply(traj2.get(i), traj1.get(j)),
                    min(FrechetHelper(i - 1, j, ca, traj1, traj2), FrechetHelper(i - 1, j - 1, ca, traj1, traj2),
                            FrechetHelper(i, j - 1, ca, traj1, traj2)));
        }
        return ca[i][j];
    }

    public enum MeasureType {
        DTW, LCSS, EDR, LORS, Hausdorff, Frechet, // actually used here
        ERP, ED, // implemented
        // TODO
        EDwP, APM, OWD, LIP, MD, STLCSS, STLC, STED, STLIP
    }
}
