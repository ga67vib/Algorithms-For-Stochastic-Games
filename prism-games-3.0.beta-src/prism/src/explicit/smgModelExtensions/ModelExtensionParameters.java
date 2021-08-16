package explicit.smgModelExtensions;

public class ModelExtensionParameters {
    public boolean extendWithBigMECs;
    public boolean extendWithProbabilityTrees;

    //BigMEC Parameters
    public int mecExtension_numMECs;
    public int mecExtension_chainLength;
    public double mecExtension_probabilityLeadingToSink;
    public double mecExtension_probabilityToGoBackToIntialState;
    public double mecExtension_probabilityOfChainSwitch;

    //ProbabilityTree Parameters
    public int probExtension_numComponents;
    public int probExtension_componentBranchingFactor;
    public int probExtension_componentTreeDepth;
    public double probExtension_probabilityLeadingToSink;
    public double probExtension_probabilityToGoBackToIntialState;
}
