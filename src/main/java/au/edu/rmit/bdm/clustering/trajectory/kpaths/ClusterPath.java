package au.edu.rmit.bdm.clustering.trajectory.kpaths;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import java.util.*;

/*
 * this store the information of each cluster
 */
public class ClusterPath {
    protected Set<Integer> clusterTrajectory; // stores all the indices of trajectories
    protected VIseries centroid; // store the centroid path
    static final int MAX_NUM_ITERATIONS = 1000; // the maximum iteration times

    // build histogram for edges and length using hashmap and Guava
    protected Multiset<Integer> edgeOcc;
    protected Multiset<Integer> lengthOcc;
    protected boolean centerChanged = true;
    int minLength = 0, maxLength = 0;
    int[] lengthAccumulation = null; // for the score computation
    Map<String, Double> checkedList = null;//store the start, end edges, and the score
    PriorityQueue<Path> queue = null;
    ArrayList<Integer> sortedFrequency = null;
    double sumEdgeOcc = 0.0;
//    int pathMinLength;

    protected int idx;// it stores the idx of the trajectory which is the centroid.
    protected double sumDistance = 0; // the sum distance in this cluster.
    protected Set<Integer> candidateList = new HashSet<>();//built from the inverted index
    protected int[] finalPath;

    public ClusterPath(int[] cl, int idx1) {
        clusterTrajectory = new HashSet<>();
        centroid = new VIseries();
        centroid.setVIseries(cl);
        if (cl != null) {
            centroid.length = cl.length;
        } else {
            centroid.length = 0;
        }
        centroid.idx = idx1;
        this.idx = idx1;
        this.finalPath = cl;
        this.edgeOcc = HashMultiset.create();
        this.lengthOcc = HashMultiset.create();
    }

    /*
     * generate the list id of trajectory that share edge with cluster
     */
    public Set<Integer> createCandidateList(Map<Integer, List<Integer>> edgeIndex, Map<Integer, int[]> dataMap) {
        if (centerChanged) {//create a new list if changed, otherwise return previous to avoid rebuild
            int[] trajectory = dataMap.get(idx);
            candidateList = new HashSet<>();
            for (int edgeId : trajectory) {
                List<Integer> trajIdList = edgeIndex.get(edgeId);
                Collections.addAll(candidateList, trajIdList.toArray(new Integer[0]));
            }
        }
        return candidateList;
    }

    /*
     * generate the list id of trajectory that share edge with cluster, or we can ignore this
     */
    public Set<Integer> createCandidateListNoDataMap(Map<Integer, List<Integer>> edgeIndex, int[] trajectory) {
        if (centerChanged) {//create a new list if changed, otherwise return previous to avoid rebuild
            candidateList = new HashSet<>();
            for (int edgeId : trajectory) {
                List<Integer> traidlist = edgeIndex.get(edgeId);
                Collections.addAll(candidateList, traidlist.toArray(new Integer[0]));
            }
        }
        return candidateList;
    }

    public Set<Integer> getCandidateList() {
        return candidateList;
    }

    /*
     * add the trajectory into the clusters
     */
    void mergeTrajectoryToCluster(ArrayList<Integer> index) {
        clusterTrajectory.addAll(index);
    }

    /*
     * add the trajectory into the clusters
     */
    void removeTrajectoryToCluster(ArrayList<Integer> index) {
        clusterTrajectory.removeAll(index);
    }

    /*
     * add the trajectory into the clusters
     */
    void addTrajectoryToCluster(int index) {
        clusterTrajectory.add(index);
    }

    /*
     * remove the trajectory into the clusters
     */
    void removeTrajectoryToCluster(int value) {
        clusterTrajectory.remove(value);
    }

    /*
     * TODO: we will try to use the data sketches instead of using hashmap
     */
    public void updateHistogramGuava(int[] tra, int idx) {
        lengthOcc.add(tra.length);
        for (int edge : tra) {
            edgeOcc.add(edge);
        }
    }

    /*
     * we will try to use the data sketches instead of using hashmap;
     */
    void removeHistogramGuava(int[] tra) {
        lengthOcc.remove(tra.length, 1);
        for (int edge : tra) {
            edgeOcc.remove(edge, 1);
        }
    }

    /*
     * we will try to use the data sketches instead of using hashmap;
     */
    public void updateHistogramGuava(Multiset<Integer> edgeH, Multiset<Integer> lengthH) {
        edgeOcc.addAll(edgeH);
        lengthOcc.addAll(lengthH);
    }

    /*
     * we will try to use the data sketches instead of using hashmap;
     */
    public void removeHistogramGuava(Multiset<Integer> edgeH, Multiset<Integer> lengthH) {
        edgeOcc.removeAll(edgeH);
        lengthOcc.removeAll(lengthH);
    }

    /*
     * return the sum distance in the cluster
     */
    public double getSumDistance() {
        return sumDistance;
    }

    /*
     * return the clusters
     */
    VIseries getClusterPath() {
        return centroid;
    }

    /*
     * return the trajectory array of this cluster
     */
    public Set<Integer> getClusterTrajectories() {
        return clusterTrajectory;
    }

    public int getTrajectoryID() {
        return idx;
    }

    public int[] getTrajectoryData() {
        return finalPath;
    }

    public boolean isCenterChanged() {
        return centerChanged;
    }

    //for bound computing
    int getFirstNSum(ArrayList<Integer> sortedFrequency, int number) {
        int sum = 0;
        for (int i = 0; i < number; i++) {
            sum += sortedFrequency.get(i);
        }
        return sum;
    }

    public void accumulateLengthOcc() {
        ArrayList<Integer> sortedTrajectoryLengths = new ArrayList<>(lengthOcc.elementSet());
        Collections.sort(sortedTrajectoryLengths);        //sorted length increasingly
        if (lengthOcc.isEmpty()) {
            System.out.println("A cluster is empty now due to bad random initialization");
            return;
        }
        minLength = sortedTrajectoryLengths.get(0);
        maxLength = sortedTrajectoryLengths.get(sortedTrajectoryLengths.size() - 1);
        lengthAccumulation = new int[maxLength - minLength + 1];
        sumEdgeOcc = 0;
        sortedFrequency = new ArrayList<>();
        for (int edge : edgeOcc.elementSet()) {//reverse way
            sortedFrequency.add(edgeOcc.count(edge));
            sumEdgeOcc += edgeOcc.count(edge);        //sorted frequency
        }
        sortedFrequency.sort(Collections.reverseOrder());
        for (int length = minLength; length <= maxLength; length++) {
            int occ = 0;
            for (int i = length; i > minLength; i--) {
                if (lengthOcc.contains(i))
                    occ += lengthOcc.count(i);
            }
            lengthAccumulation[length - minLength] = occ;
        }
    }

    /*
     * We use Multiset in the google library to update the histogram in a faster way by scanning each trajectory
     */
    double extractNewPathGuava(Map<Integer, int[]> dataMap, RunLog runRecord,
                               Map<Integer, Integer> traLengthMap, Map<Integer, Integer> trajectoryHistogram) {
        int min = Integer.MAX_VALUE;
        int trajId = 0;
        accumulateLengthOcc();
        for (int clusterTrajIdx : clusterTrajectory) {
            // read each trajectory in this cluster,
            // it is slow as we need to read every trajectory, we can construct a weighted graph and choose the route
            int traLength = traLengthMap.get(clusterTrajIdx);
            int bound = Math.min(trajectoryHistogram.get(clusterTrajIdx), getFirstNSum(sortedFrequency, traLength));
            int sumLength = edgeOcc.size();
            int i = minLength;
            while (i <= traLength) {
                sumLength += lengthAccumulation[i - minLength];
                i++;
            }
            if (sumLength - bound >= min) {
                continue;// skip reading the trajectory data
            }
            int[] tra = dataMap.get(clusterTrajIdx);
            int sumFrequency = 0;
            for (int edge : tra) {        //compute the vertex frequency of each trajectory
                sumFrequency += edgeOcc.count(edge);
            }
            int sumDis = sumLength - sumFrequency;
            if (min > sumDis) {//choose the one which can minimize the sum of distance
                min = sumDis;
                trajId = clusterTrajIdx;
            }
        }
        if (idx == trajId) {
            centerChanged = false;//center does not change;
        }
        sumDistance = min;
        int[] a = dataMap.get(trajId);
        //	System.out.println(edgeOcc.size()+" "+min+" "+a.length);
        double drift = 0;
        if (centerChanged) {
            drift = Intersection(finalPath, a, finalPath.length, a.length);
        }
        centroid.setVIseries(a);
        centroid.idx = trajId;
        idx = trajId;
        finalPath = a;
        return drift;
    }

    /* A greedy algorithm when the optimal result does not change any more and better than last result,
     * we will stop exploring when the length is close to the maximum length it can be
     */
    public double extractNewPathFrequency(HashMap<Integer, ArrayList<Integer>> forwardGraph,
                                          HashMap<Integer, ArrayList<Integer>> backwardGraph) {
        checkedList = new HashMap<>();
        queue = new PriorityQueue<>();
        ArrayList<Integer> optimal = new ArrayList<>();
        accumulateLengthOcc();
        double optimalScore = Double.MAX_VALUE;

        if (finalPath != null) {
            centerChanged = false;
            optimalScore = computeScore(finalPath);//initialized to min score using last centroid
            for (int a : finalPath)
                optimal.add(a);
            double lowerBound = estimateLowerBound(optimalScore, optimal.size());
            if (optimalScore > lowerBound) {
                Path aPath = new Path(optimal, optimalScore, lowerBound);
                queue.add(aPath);
            }
        }

        for (int edge : edgeOcc.elementSet()) {//initialize the edges as paths
            ArrayList<Integer> arrayList = new ArrayList<>();
            arrayList.add(edge);
            double score = computeScore(arrayList);
            double lowerBound = estimateLowerBound(score, 1);
            if (optimalScore > lowerBound) {
                Path aPath = new Path(arrayList, score, lowerBound);
                queue.add(aPath);
            }
        }

        int iter = 0;
        while (!queue.isEmpty()) {
            iter++;
            Path candidate = queue.poll();
            double lowerBound = candidate.getLowerbound();// compute the score
            double score = candidate.getScore();
            if (lowerBound > optimalScore || iter >= MAX_NUM_ITERATIONS) {//termination as all possible has been checked
                break;
            }
            ArrayList<Integer> Can = candidate.getPath();
            int start = Can.get(0);
            int end = Can.get(Can.size() - 1);// the last edge
            if (backwardGraph.containsKey(start)) {
                ArrayList<Integer> backAppend = backwardGraph.get(start);
                if (backAppend != null)
                    for (int ids : backAppend) {
                        if (edgeOcc.contains(ids) && !Can.contains(ids)) {
                            ArrayList<Integer> newCan = new ArrayList<>(Can);
                            newCan.add(0, ids);    //insert to the beginning
                            if (forwardGraph.containsKey(end)) {
                                ArrayList<Integer> startAppend1 = forwardGraph.get(end);// no circle
                                if (startAppend1 != null && startAppend1.contains(ids)) {
                                    continue;
                                }
                            }
                            double newscore = computeScoreWithPrevious(ids, score, newCan.size());
                            if (newscore < optimalScore) {
                                centerChanged = true;//the center has changed
                                optimalScore = newscore;
                                optimal = newCan;
                            }
                            lowerBound = estimateLowerBound(newscore, newCan.size());
                            String signature = signaturePath(ids, end);
                            if (checkWithLowerBound(signature, lowerBound)) {
                                continue;
                            }
                            if (lowerBound < optimalScore) {
                                queue.add(new Path(newCan, newscore, lowerBound));
                            }
                        }
                    }
            }

            if (forwardGraph.containsKey(end)) {
                ArrayList<Integer> startAppend = forwardGraph.get(end);
                if (startAppend != null)
                    for (int ids : startAppend) {
                        if (edgeOcc.contains(ids) && !Can.contains(ids)) {//connected and no repetitive edges
                            ArrayList<Integer> newCan = new ArrayList<>(Can);
                            newCan.add(ids);
                            if (forwardGraph.containsKey(ids)) {
                                ArrayList<Integer> startAppend1 = forwardGraph.get(ids);
                                if (startAppend1 != null && startAppend1.contains(start)) {
                                    continue;
                                }
                            }
                            double newScore = computeScoreWithPrevious(ids, score, newCan.size());
                            if (newScore < optimalScore) {
                                centerChanged = true;
                                optimalScore = newScore;
                                optimal = newCan;
                            }
                            lowerBound = estimateLowerBound(newScore, newCan.size());
                            String signature = signaturePath(start, ids);
                            if (checkWithLowerBound(signature, lowerBound))    //repetitive path
                                continue;
                            if (lowerBound < optimalScore) {//if this path is promising
                                queue.add(new Path(newCan, newScore, lowerBound));
                            }
                        }
                    }
            }
        }
        sumDistance = optimalScore;
        //	System.out.println(optimalScore+" "+optimal.size());
        Collections.sort(optimal);
        int[] centroidData = optimal.stream().mapToInt(i -> i).toArray();
        double drift = 0;
        if (finalPath != null && centerChanged)
            drift = Intersection(finalPath, centroidData, finalPath.length, centroidData.length);
        centroid.setVIseries(centroidData);
        finalPath = centroidData;
        return drift;
    }

    /*
     * the data needs to be sorted before the intersection
     */
    public int Intersection(int[] arr1, int[] arr2, int m, int n) {
        int i = 0, j = 0;
        int dist = 0;
        while (i < m && j < n) {
            if (arr1[i] < arr2[j])
                i++;
            else if (arr2[j] < arr1[i])
                j++;
            else {
                dist++;
                i++;
                j++;
            }
        }
        return Math.max(m, n) - dist;
    }

    /*
     * check whether it exists in the checked list
     */
    public boolean checkWithLowerBound(String aString, double lowerBound) { // TODO ?
        if (checkedList.containsKey(aString)) {
            double previous = checkedList.get(aString);
            if (lowerBound >= previous) {
                return false;
            } else {
                checkedList.put(aString, lowerBound);// add to the list
                return true;
            }
        } else {
            checkedList.put(aString, lowerBound);// add to the list
            return true;
        }
    }

    /*
     * compute the objective score based on the edge histogram and length histogram,
     */
    public double computeScore(int[] path) {
        double weight = sumEdgeOcc;
        for (int value : path) {
            weight -= edgeOcc.count(value);
        }
        if (path.length > minLength) {
            for (int i = 1; i <= Math.min(path.length - minLength, maxLength - minLength); i++)
                weight += lengthAccumulation[i];
        }
        return weight;
    }

    /*
     * compute the objective score based on the edge histogram and length histogram
     */
    public double computeScore(ArrayList<Integer> path) {
        return computeScore(path.stream().mapToInt(i -> i).toArray());
    }

    /*
     * compute the objective score based on the edge histogram and length histogram
     */
    public double computeScoreWithPrevious(int newEdge, double weight, int length) {
        weight -= edgeOcc.count(newEdge);
        if (length > minLength) {
            weight += lengthAccumulation[length - minLength];
        }
        return weight;
    }

    /*
     * convert to a short unique string using the start and end edge, and the length,
     * this can be used as a dominance table,
     */
    public String signaturePath(int a, int b) {
        return a + "_" + b;
    }

    //add the rest heavy edges to compute the lower bound of the score
    public double estimateLowerBound(double score, int length) {
        double lowerBound = score;
        int i = length;
        while (i < minLength && i - length < sortedFrequency.size()) {
            lowerBound -= sortedFrequency.get(i - length);
            i++;
        }
        while (i <= maxLength && i - length < sortedFrequency.size()
                && lengthAccumulation[i - minLength] < sortedFrequency.get(i - length)) {
            lowerBound += lengthAccumulation[i - minLength] - sortedFrequency.get(i - length);
            i++;
        }
        return lowerBound;
    }

    /*
     * the maximum length of the final output path
     */
    public int estimateMax() {
        int i = minLength;
        while (lengthAccumulation[i - minLength] < sortedFrequency.get(i - minLength)) {
            i++;
        }
        return i;
    }
}
