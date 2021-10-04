package explicit.smgModelExtensions;

import explicit.STPGExplicit;

public class ModelExtensionResult {
    public int addedInitialState;
    public STPGExplicit resultingSTPG;

    public ModelExtensionResult(int addedInitialState, STPGExplicit resultingSTPG) {
        this.addedInitialState = addedInitialState;
        this.resultingSTPG = resultingSTPG;
    }
}