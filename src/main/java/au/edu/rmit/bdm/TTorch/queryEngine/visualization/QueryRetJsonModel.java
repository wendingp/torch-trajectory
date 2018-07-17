package au.edu.rmit.bdm.TTorch.queryEngine.visualization;

import au.edu.rmit.bdm.TTorch.base.Torch;
import au.edu.rmit.bdm.TTorch.queryEngine.query.QueryResult;

import java.util.ArrayList;
import java.util.List;

class QueryRetJsonModel {

    final String queryType;
    final boolean mappingSucceed;
    final TrajJsonModel raw;
    final TrajJsonModel mapped;
    final int retSize;
    final List<TrajJsonModel> ret;

    QueryRetJsonModel(QueryResult queryResult){
        this.queryType = queryResult.queryType;
        this.mappingSucceed = queryResult.mappingSucceed;
        this.raw = queryType.equals(Torch.QueryType.RangeQ) ? null : new TrajJsonModel(queryResult.rawQuery);
        this.mapped = queryType.equals(Torch.QueryType.RangeQ) ? null : new TrajJsonModel(queryResult.mappedQuery);
        this.retSize = queryResult.ret.size();
        this.ret = new ArrayList<>(queryResult.ret.size());
        ret.addAll(Formater.model(queryResult.ret));
    }

}