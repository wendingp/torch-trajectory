package au.edu.rmit.bdm.Torch.queryEngine.similarity;

import au.edu.rmit.bdm.Torch.base.helper.GeoUtil;
import au.edu.rmit.bdm.Torch.base.model.TrajEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

import static au.edu.rmit.bdm.Torch.base.helper.GeoUtil.min;


/**
 * Distances
 * usage:
 * DistanceFunction<TorVertex, TorVertex> distFunc = (p1, p2) -> GeoUtil.distance(p1, p2);
 * Comparator<TorVertex> comparator = (p1, p2) -> {
 * double dist = GeoUtil.distance(p1, p2);
 * if (dist < 8) return 0;
 * return 1;
 * };
 * SimilarityFunction<TorVertex> SIM_MEASURE = new SimilarityFunction<>(distFunc, comparator);
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
            if (GeoUtil.distance(p1, p2) <= EPSILON) return 0;
            return 1;
        };

        DEFAULT = new SimilarityFunction<>(distFunc, comparator);
    }

    private final DistanceFunction distFunc;

    public Comparator<T> comparator;
//    public String measure;

    public SimilarityFunction(DistanceFunction<T, T> distFunc, Comparator<T> comparator) {
        this.distFunc = distFunc;
        this.comparator = comparator;
    }

    /**
     * ED
     *
     * @param t1
     * @param t2
     * @return
     */
    public double EuclideanDistance(List<T> t1, List<T> t2) { // ED O(mn)
        if (t1.size() > t2.size()) {
            List<T> tmp = t1;
            t1 = t2;
            t2 = tmp;
        }
        int n = t1.size(), m = t2.size();
        double minDistance = Double.MAX_VALUE;
        for (int st = 0; st + n <= m; ++st) {
            minDistance = Math.min(minDistance, EuclideanDistanceSameLength(t1, t2.subList(st, st + n)));
        }
        return minDistance;
    }

    private double EuclideanDistanceSameLength(List<T> t1, List<T> t2) { // ED of 2 same size trajectories O(n)
        if (t1.size() != t2.size()) {
            throw new IllegalArgumentException("T1 should be of the same length as T2");
        }
        double dist = 0.0;
        for (int i = 0; i < t1.size(); ++i) {
            dist += Math.sqrt(distFunc.apply(t1.get(i), t2.get(i)));
        }
        return dist / t1.size(); // TODO should be average here?
//        return dist;
    }

    public double LongestCommonSubsequence(List<T> t1, List<T> t2, int theta) { // LCSS
        if (t1 == null || t2 == null || t1.size() == 0 || t2.size() == 0)
            return 0;

        int[][] dp = new int[t1.size()][t2.size()];

        if (comparator.compare(t1.get(0), t2.get(0)) == 0)
            dp[0][0] = 1;

        for (int i = 1; i < t1.size(); ++i) {
            if (comparator.compare(t1.get(i), t2.get(0)) == 0) {
                dp[i][0] = 1;
            } else {
                dp[i][0] = dp[i - 1][0];
            }
        }

        for (int i = 1; i < t2.size(); ++i) {
            if (comparator.compare(t2.get(i), t1.get(0)) == 0) {
                dp[0][i] = 1;
            } else {
                dp[0][i] = dp[0][i - 1];
            }
        }

        for (int i = 1; i < t1.size(); ++i) {
            for (int j = 1; j < t2.size(); ++j) {
                if (Math.abs(i - j) <= theta) { // TODO ?
                    if (comparator.compare(t1.get(i), t2.get(j)) == 0) {
                        dp[i][j] = 1 + dp[i - 1][j - 1];
                    } else {
                        dp[i][j] = Math.max(dp[i - 1][j], dp[i][j - 1]);
                    }
                }
            }
        }
        return dp[t1.size() - 1][t2.size() - 1];
    }

    /**
     * ERP
     *
     * @param t1
     * @param t2
     * @param g
     * @return
     */
    public double EditDistanceWithRealPenalty(List<T> t1, List<T> t2, T g) { // TODO what is g?
        if (t1 == null || t1.size() == 0) { // TODO why allow null?
            double res = 0.0;
            if (t2 != null) {
                for (T t : t2) {
                    res += distFunc.apply(t, g);
                }
            }
            return res;
        }
        if (t2 == null || t2.size() == 0) {
            double res = 0.0;
            for (T t : t1) {
                res += distFunc.apply(t, g);
            }
            return res;
        }

        double[][] dp = new double[t1.size() + 1][t2.size() + 1];
        for (int i = 1; i <= t1.size(); ++i) {
            dp[i][0] = distFunc.apply(t1.get(i - 1), g) + dp[i - 1][0];
        }
        for (int j = 1; j <= t2.size(); ++j) {
            dp[0][j] = distFunc.apply(t2.get(j - 1), g) + dp[0][j - 1];
        }

        for (int i = 1; i <= t1.size(); ++i) {
            for (int j = 1; j <= t2.size(); ++j) {
                dp[i][j] = min(dp[i - 1][j - 1] + distFunc.apply(t1.get(i - 1), t2.get(j - 1)),
                        dp[i - 1][j] + distFunc.apply(t1.get(i - 1), g),
                        dp[i][j - 1] + distFunc.apply(g, t2.get(j - 1)));
            }
        }
        return dp[t1.size()][t2.size()] / Math.max(t1.size(), t2.size());
    }

    /**
     * EDR
     *
     * @param t1
     * @param t2
     * @return
     */
    public double EditDistanceOnRealSequence(List<T> t1, List<T> t2) {
        if (t1.isEmpty() || t2.isEmpty()) {
            if (t1.isEmpty() && t2.isEmpty()) return 0;
            return Math.max(t1.size(), t2.size());
        }
        int[][] dp = new int[t1.size() + 1][t2.size() + 1];

        for (int i = 1; i <= t1.size(); ++i) {
            dp[i][0] = i;
        }
        for (int j = 1; j <= t2.size(); ++j) {
            dp[0][j] = j;
        }

        for (int i = 1; i <= t1.size(); ++i) {
            for (int j = 1; j <= t2.size(); ++j) {
                int subCost = 1;
                if (comparator.compare(t1.get(i - 1), t2.get(j - 1)) == 0) {
                    subCost = 0;
                }
                dp[i][j] = min(dp[i - 1][j - 1] + subCost, dp[i - 1][j] + 1, dp[i][j - 1] + 1);
            }
        }
        return dp[t1.size()][t2.size()];
    }

    /**
     * DTW
     *
     * @param t1
     * @param t2
     * @return
     */
    public double DynamicTimeWarping(List<T> t1, List<T> t2) {
        if (t1.size() == 0 && t2.size() == 0) return 0;
        if (t1.size() == 0 || t2.size() == 0) return Integer.MAX_VALUE;

        double[][] dp = new double[t1.size() + 1][t2.size() + 1];

        for (int i = 1; i <= t1.size(); ++i) {
            dp[i][0] = Integer.MAX_VALUE;
        }
        for (int j = 1; j <= t2.size(); ++j) {
            dp[0][j] = Integer.MAX_VALUE;
        }

        for (int i = 1; i <= t1.size(); ++i) {
            for (int j = 1; j <= t2.size(); ++j) {
                dp[i][j] = distFunc.apply(t1.get(i - 1), t2.get(j - 1)) + min(dp[i - 1][j - 1], dp[i - 1][j], dp[i][j - 1]);
            }
        }
        return dp[t1.size()][t2.size()];
    }

    /**
     * Hausdorff
     *
     * @param t1
     * @param t2
     * @return
     */
    public double Hausdorff(List<T> t1, List<T> t2) { // refactored to O(mn)
        double[][] distMatrix = new double[t2.size()][t1.size()]; // TODO why t2 then t1 ?
        for (int i = 0; i < distMatrix.length; ++i) {
            for (int j = 0; j < distMatrix[0].length; ++j) {
                distMatrix[i][j] = distFunc.apply(t1.get(j), t2.get(i));
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
     * @param t1
     * @param t2
     * @return
     */
    public double Frechet(List<T> t1, List<T> t2) {
        double[][] ca = new double[t2.size()][t1.size()]; // TODO ?
        for (int i = 0; i < t2.size(); ++i) {
            for (int j = 0; j < t1.size(); ++j) {
                ca[i][j] = -1.0D;
            }
        }
        return FrechetHelper(t2.size() - 1, t1.size() - 1, ca, t1, t2);
    }

    private double FrechetHelper(int i, int j, double[][] ca, List<T> t1, List<T> t2) {
        if (ca[i][j] > -1.0D)
            return ca[i][j];
        if (i == 0 && j == 0) {
            ca[i][j] = distFunc.apply(t1.get(0), t2.get(0));
        } else if (j == 0) {
            ca[i][j] = Math.max(FrechetHelper(i - 1, 0, ca, t1, t2), distFunc.apply(t2.get(i), t1.get(0)));
        } else if (i == 0) {
            ca[i][j] = Math.max(FrechetHelper(0, j - 1, ca, t1, t2), distFunc.apply(t2.get(0), t1.get(j)));
        } else {
            ca[i][j] = Math.max(distFunc.apply(t2.get(i), t1.get(j)), min(
                    FrechetHelper(i - 1, j, ca, t1, t2),
                    FrechetHelper(i - 1, j - 1, ca, t1, t2),
                    FrechetHelper(i, j - 1, ca, t1, t2)
            ));
        }
        return ca[i][j];
    }

    public enum MeasureType {
        DTW, LCSS, EDR, LORS, Hausdorff, Frechet, // actually used here
        ERP, ED, // implemented
        // TODO
        EDwP, APM, OWD, LIP,
        MD, STLCSS, STLC, STED, STLIP
    }
}
