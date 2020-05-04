package au.edu.rmit.bdm.clustering.trajectory.streaming;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

/*
 * it stores the inverted index for each edge
 */
public class streamEdge {
    TreeMap<Integer, ArrayList<Integer>> timeIndex;//sorted timestamp, the carid list
    int id;//the id of each edge

    public streamEdge(int id) {
        this.id = id;
        timeIndex = new TreeMap<>();
    }

    public void addCars(int cartime, int[] carsid) {
        ArrayList<Integer> traidlist;
        if (timeIndex.containsKey(cartime)) {
            traidlist = timeIndex.get(cartime);
        } else {
            traidlist = new ArrayList<>();
        }
        for (int i : carsid) {
            traidlist.add(i);
        }
        timeIndex.put(cartime, traidlist);
    }

    /*
     * delete expired data which will not be used for data and index within in the sliding window
     */
    public void removeExpired(int expiredtime, int formerexpired, Map<Integer, ArrayList<Integer>> dataset_remove) {
        if (timeIndex.isEmpty()) return;
        for (int timeId : timeIndex.keySet()) {
            if (timeId < expiredtime && timeId >= formerexpired) {
                ArrayList<Integer> cars = timeIndex.get(timeId);
                for (int carId : cars) {
                    ArrayList<Integer> edgeIds;
                    if (dataset_remove.containsKey(carId)) {
                        edgeIds = dataset_remove.get(carId);
                    } else {
                        edgeIds = new ArrayList<>();
                    }
                    edgeIds.add(id);
                    dataset_remove.put(carId, edgeIds);
                }
            }
        }
    }
}
