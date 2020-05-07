package au.edu.rmit.bdm.clustering.trajectory.kpaths;

/*
 * this class is used for recording the running statistics
 */
public class RunLog {
    private double assignmentTime;
    private double refinementTime;
    private int iterationTimes;
    private double simiComputationTime;
    private double ioTime;
    private double histogramTime;
    private double accumulatedTime;

    public RunLog() {
        assignmentTime = 0;
        refinementTime = 0;
        iterationTimes = 0;
        simiComputationTime = 0;
        ioTime = 0;
        histogramTime = 0;
        accumulatedTime = 0;
    }

    public void clear() {
        assignmentTime = 0;
        refinementTime = 0;
        iterationTimes = 0;
        simiComputationTime = 0;
        ioTime = 0;
        histogramTime = 0;
        accumulatedTime = 0;
    }

    public void addAssignmentTime(double time) {
        assignmentTime += time;
    }

    public void addRefinementTime(double time) {
        refinementTime += time;
    }

    public void addsimiComputationTime(double time) {
        simiComputationTime += time;
    }

    public void addIOTime(double time) {
        ioTime += time;
    }

    public void addHistogramTime(double time) {
        histogramTime += time;
    }

    public void addAccumulatedTime(double time) {
        accumulatedTime += time;
    }

    public void setIterationtimes(int ite) {
        iterationTimes = ite;
    }

    public double getAssignmentTime() {
        return assignmentTime;
    }

    public double getRefinementTime() {
        return refinementTime;
    }

    public double getsimiComputationTime() {
        return simiComputationTime;
    }

    public double getAlltime() {
        return assignmentTime + refinementTime;
    }

    public void printLog() {
        double alltime = assignmentTime + refinementTime;
        System.out.println("\n==============================================================");
        System.out.println("it used " + iterationTimes + " iterations to find the centroids");
        System.out.println("the overall running time: " + alltime + "s");
        System.out.println("\tthe assignment time: " + (assignmentTime - histogramTime) + "s");
        System.out.println("\t\tthe IO and Index access time: " + ioTime + "s");
        System.out.println("\t\tthe similarity computation time: " + simiComputationTime + "s");
        System.out.println("\tthe refinement time: " + (refinementTime + histogramTime) + "s");
        System.out.println("\t\tthe edge histogram building time: " + histogramTime + "s");
        System.out.println("\t\tthe accumulated histogram time: " + accumulatedTime + "s");
        System.out.println("==============================================================\n");
    }
}
