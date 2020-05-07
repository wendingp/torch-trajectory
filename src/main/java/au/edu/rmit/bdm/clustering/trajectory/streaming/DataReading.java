package au.edu.rmit.bdm.clustering.trajectory.streaming;

import au.edu.rmit.bdm.clustering.trajectory.kpaths.Util;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;

public class DataReading extends Thread {

    private static Map<String, Integer> edgeMapping; // the key stores the composition of camera id, lane id, direction, the value stores the new id
    private static Map<String, Integer> vehicleMapping;// the key stores the original plate number, the value is the new id.
    // static int globalStartTime = 151792;


    /*
     * convert the raw data from oracle to a simple version
     */
    public static void convertToEdges(String path, String output, String newEdgeFile, String newCarFile) {
        int edgeCounter = 1;
        int carCounter = 1;
        edgeMapping = new HashMap<>();
        vehicleMapping = new HashMap<>();
        try {
            Scanner in = new Scanner(new BufferedReader(new FileReader(path)));
            while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
                String str = in.nextLine();
                String strTrim = str.trim();
                String[] record = strTrim.split(",");
                if (record[0].equals("LKBH") || record.length < 6)
                    continue; //|| StringUtils.isNumeric(record[5])==false
                String oldEdge = record[0] + "_" + record[1] + "_" + record[2];//can be combined with the lane number for a high granularity
                int newEdge;
                if (edgeMapping.containsKey(oldEdge)) {
                    newEdge = edgeMapping.get(oldEdge);
                } else {
                    newEdge = edgeCounter++;
                    edgeMapping.put(oldEdge, newEdge);
                }
                int newCar;
                if (vehicleMapping.containsKey(record[3])) {
                    newCar = vehicleMapping.get(record[3]);
                } else {
                    newCar = carCounter++;
                    vehicleMapping.put(record[3], newCar);
                }
                double normalizedTime = Double.parseDouble(record[5]);
                String newRecord = (int) normalizedTime + "," + newCar + "," + newEdge + "\n";
                System.out.println(newRecord);
                Util.write(output, newRecord);
            }
            in.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        for (String oldCar : vehicleMapping.keySet()) {
            int newCar = vehicleMapping.get(oldCar);
            String newRecord = newCar + "," + oldCar + "\n";
            Util.write(newCarFile, newRecord);
        }
        for (String oldEdge : edgeMapping.keySet()) {
            int newEdge = edgeMapping.get(oldEdge);
            String newRecord = newEdge + "," + oldEdge + "\n";
            Util.write(newEdgeFile, newRecord);
        }
    }

    /*
     * convert to the standard format, each line is the timestamp; edge: carid1, carid2...;
     */
    public static void convertStandardFormat(String path, String output) {
        Map<Integer, Set<Integer>> storeSecondRecord = null;
        int tempid = Integer.MIN_VALUE;
        try {
            Scanner in = new Scanner(new BufferedReader(new FileReader(path)));
            while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
                String str = in.nextLine();
                String strr = str.trim();
                String[] record = strr.split(",");
                int cartime = Integer.parseInt(record[0]);
                int carid = Integer.parseInt(record[1]);
                int edgeid = Integer.parseInt(record[2]);
                if (tempid != cartime) {
                    if (storeSecondRecord != null) {
                        int time = tempid + 12;
                        StringBuilder content = new StringBuilder(time + ";");
                        for (int edgeid1 : storeSecondRecord.keySet()) {
                            content.append(edgeid1).append(":");
                            Set<Integer> tralist = storeSecondRecord.get(edgeid1);
                            for (int carid1 : tralist)
                                content.append(carid1).append(",");
                            content = new StringBuilder(content.substring(0, content.length() - 1));
                            content.append(";");
                        }
                        content = new StringBuilder(content.substring(0, content.length() - 1));
                        Util.write(output, content + "\n");
                    }
                    storeSecondRecord = new TreeMap<>();
                    tempid = cartime;
                }
                Set<Integer> tralist;
                assert storeSecondRecord != null;
                if (storeSecondRecord.containsKey(edgeid)) {
                    tralist = storeSecondRecord.get(edgeid);
                } else {
                    tralist = new HashSet<>();
                }
                tralist.add(carid);
                storeSecondRecord.put(edgeid, tralist);
            }
            in.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
