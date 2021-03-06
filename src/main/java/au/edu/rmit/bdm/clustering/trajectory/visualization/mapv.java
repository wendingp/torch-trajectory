package au.edu.rmit.bdm.clustering.trajectory.visualization;

import java.io.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.Scanner;

public class mapv {
    private static Scanner in;

    /*
     * show the specific path in the road
     */
    public static void generateClusterPath(Map<Integer, int[]> map, Map<Integer, String> edgeInfo,
                                           ArrayList<Integer> clusterIDs, String output) {
        write(output, "geometry\n");
        for (Integer clusterID : clusterIDs) {
            StringBuilder content = new StringBuilder("\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [");
            int[] cluster = map.get(clusterID);
            // System.out.println(Arrays.toString(cluster));
            for (int value : cluster) {
                String singleEdge = edgeInfo.get(value);
                String[] abc = singleEdge.split(",");
                for (int i = 0; i < abc.length / 2; i++) {
                    content.append("[").append(abc[i + abc.length / 2]).append(",").append(abc[i]).append("],");
                }
            }
            content = new StringBuilder(content.substring(0, content.length() - 1));
            content.append("]}\"");
            write(output, content + "\n");
        }
    }

    /*
     * show the specific path in the road
     */
    public static void generateClusterPath1(Map<Integer, int[]> map, Map<Integer, String> edgeInfo,
                                            ArrayList<int[]> clusterIDs, String output) {
        write(output, "geometry\n");
        for (int[] cluster : clusterIDs) {
            for (int value : cluster) {
                String single_edge = edgeInfo.get(value);
                StringBuilder content = new StringBuilder("\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [");
                String[] abc = single_edge.split(",");
                for (int i = 0; i < abc.length / 2; i++) {
                    content.append("[").append(abc[i + abc.length / 2]).append(",").append(abc[i]).append("],");
                }
                content = new StringBuilder(content.substring(0, content.length() - 1));
                content.append("]}\"");
                write(output, content + "\n");
            }
        }
    }

    /*
     * show the specific sorted path in the road network
     */
    public static void generateClusterPathSorted(Map<Integer, int[]> map, Map<Integer, String> edgeInfo,
                                                 ArrayList<Integer> clusterIDs, String output) {
        write(output, "geometry\n");
        for (Integer clusterID : clusterIDs) {
            int[] cluster = map.get(clusterID);
            // System.out.println(Arrays.toString(cluster));
            for (int value : cluster) {
                StringBuilder content = new StringBuilder("\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [");
                String single_edge = edgeInfo.get(value);
                String[] abc = single_edge.split(",");
                for (int i = 0; i < abc.length / 2; i++) {
                    content.append("[").append(abc[i + abc.length / 2]).append(",").append(abc[i]).append("],");
                }
                content = new StringBuilder(content.substring(0, content.length() - 1));
                content.append("]}\"");
                write(output, content + "\n");
            }
        }
    }

    /*
     * show the subset of edges
     */
    public static void generateHighEdges(Map<Integer, String> edgeInfo,
                                         ArrayList<Integer> edgeIDs, String output) {
        write(output, "geometry\n");
        for (Integer edgeID : edgeIDs) {
            StringBuilder content = new StringBuilder("\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [");
            String single_edge = edgeInfo.get(edgeID);
            String[] abc = single_edge.split(",");
            for (int i = 0; i < abc.length / 2; i++) {
                content.append("[").append(abc[i + abc.length / 2]).append(",").append(abc[i]).append("],");
            }
            content = new StringBuilder(content.substring(0, content.length() - 1));
            content.append("]}\"");
            write(output, content + "\n");
        }
    }


    /*
     * generate the graph for mapv based on the edge graph data
     */
    public static void generate_mapv_graph_porto(String edge, String output) {
        write(output, "geometry\n");
        try {
            in = new Scanner(new BufferedReader(new FileReader(edge)));
            while (in.hasNextLine()) {
                StringBuilder content = new StringBuilder("\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [");
                String str = in.nextLine();
                String strr = str.trim();
                String[] abc = strr.split(",");
                for (int i = 0; i < abc.length / 2; i++) {
                    content.append("[").append(abc[i + abc.length / 2]).append(",").append(abc[i]).append("]");
                    if (i < abc.length / 2 - 1)
                        content.append(",");
                }
                content.append("]}\"");
                write(output, content + "\n");
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    // write the information into files
    public static void write(String fileName, String content) {
        RandomAccessFile randomFile = null;
        try {
            randomFile = new RandomAccessFile(fileName, "rw");
            long fileLength = randomFile.length();
            randomFile.seek(fileLength);
            randomFile.writeBytes(content);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (randomFile != null) {
                try {
                    randomFile.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /*
     * convert the sorted sequence into a path in the road network.
     */
    public static void reconvertTrajectory() {
        throw new UnsupportedOperationException();
    }
}