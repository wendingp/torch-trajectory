package au.edu.rmit.bdm.Torch.queryEngine.visualization;

import au.edu.rmit.bdm.Torch.base.model.TrajEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

class Geometry {
    private static final Logger logger = LoggerFactory.getLogger(Geometry.class);

    // String type = "LineString";

    //lng, lat
    //order is vital
    double[][] coordinates;

    Geometry(List<TrajEntry> path) {
        if (path == null) {
            logger.warn("path for Geometry is null");
            return;
        }
        coordinates = new double[path.size()][2];
        for (int i = 0; i < path.size(); i++) {
            coordinates[i][0] = path.get(i).getLng();
            coordinates[i][1] = path.get(i).getLat();
        }
    }
}
