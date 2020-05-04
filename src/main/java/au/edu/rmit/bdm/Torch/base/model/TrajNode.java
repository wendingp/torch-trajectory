package au.edu.rmit.bdm.Torch.base.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A trajNode models a node on a trajectory.
 * Besides GPS coordinate, each node could also carries other information such as timestamp, current speed etc.
 * However, algorithm from raw trajectory nodes to graph vertices for other kind of information requires further
 * research and experiment.
 * This function is expected to be developed in future.
 */
public class TrajNode extends Coordinate {
    private final long INVALID_TIME_VAL = -1;
    public int id;
    private long _time = INVALID_TIME_VAL;
    List<String> _bundle = new ArrayList<>(3);

    public TrajNode(double lat, double lng) {
        super(lat, lng);
    }

    public TrajNode(double lat, double lng, long time) {
        this(lat, lng);
        _time = time;
    }

    public void addExtraInfo(String extra) {
        _bundle.add(extra);
    }

    public void setTime(long time) {
        _time = time;
    }

    public long getTime() {
        if (_time == INVALID_TIME_VAL) {
            throw new IllegalStateException("try to fetch time in node not containing time information");
        }
        return _time;
    }

    @Override
    public int getId() {
        return id;
    }
}
