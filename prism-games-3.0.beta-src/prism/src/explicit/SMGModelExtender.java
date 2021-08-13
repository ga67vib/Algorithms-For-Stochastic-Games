package explicit;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import parser.State;
import prism.PrismException;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;

public class SMGModelExtender {
    public static STPG extendModel(STPG stpg, BitSet remain, String configFilePath) throws PrismException {

        ModelExtensionParameters config = ConfigParserJSON.read(configFilePath);

        Iterable<Integer> oldInitialStatesIterable = stpg.getInitialStates();
        ArrayList<Integer> oldInitialStates = new ArrayList<>();
        for (Integer initialState : oldInitialStatesIterable) {
            oldInitialStates.add(initialState);
        }

        STPGExplicit newStpg = (STPGExplicit) stpg;//(STPGExplicit) createCopyOfExistingSTPG(stpg);
        ArrayList<Integer> newInitialStates = new ArrayList<>();

        int sink = newStpg.addState();
        Distribution sinkSelfLoop = new Distribution();
        sinkSelfLoop.add(sink, 1.0);
        newStpg.addChoice(sink, sinkSelfLoop);
        if (remain != null) remain.set(sink);
        addToStateList(newStpg, sink, sink);

        Iterable<Integer> initialStates;

        if (config.extendWithMECs) {
            // Iterate over every possible startState for BigMEC
            initialStates = newStpg.getInitialStates();
            for (int initialState : initialStates) {
                ModelExtensionResult modelExtensionResult = addStartingMEC(newStpg, remain, initialState, config.mecExtension_numMECs, config.mecExtension_chainLength, config.mecExtension_probabilityLeadingToSink, sink, config.mecExtension_probabilityOfChainSwitch);
                newStpg = modelExtensionResult.resultingSTPG;
                newInitialStates.add(modelExtensionResult.addedInitialState);
            }

            newStpg.clearInitialStates();
            for (int newInitialState : newInitialStates) {
                newStpg.addInitialState(newInitialState);
            }
        }

        if (config.extendWithProbabilityComponents) {
            // Iterate over every possible startState for BigMEC
            initialStates = newStpg.getInitialStates();
            for (int initialState : initialStates) {
                ModelExtensionResult modelExtensionResult = addStartingProbabilityGraph(newStpg, remain, initialState, config.probExtension_numComponents, config.probExtension_componentTreeDepth, config.probExtension_componentBranchingFactor,
                        config.probExtension_probabilityToGoBackToIntialState, config.probExtension_probabilityLeadingToSink, sink);
                newStpg = modelExtensionResult.resultingSTPG;
                newInitialStates.add(modelExtensionResult.addedInitialState);
            }

            newStpg.clearInitialStates();
            for (int newInitialState : newInitialStates) {
                newStpg.addInitialState(newInitialState);
            }
        }
        return newStpg;
    }

    /**
     * Adds MECs in front of the initial states. The MEC is of type "BigMEC", since one may not want to have as many
     * probabilistic edges as in ManyMECs.
     *
     * @param stpg
     * @param initialState:         Before which state should MECs be chained
     * @param numChainedMECs:       How many MECs are there before initial State?
     * @param chainLengthInMEC:     How long should one MEC be? (Total size of MEC is 2 * chainLengthInMEC + 1)
     * @param sinkProbabilityPerMEC: Probability of reaching a sink on the actions that leave the component
     * @param probabilityInMEC:     Probability of moving to another chain
     */
    private static ModelExtensionResult addStartingMEC(STPGExplicit stpg, BitSet remain, int initialState, int numChainedMECs, int chainLengthInMEC,
                                                       double sinkProbabilityPerMEC, int sinkState, double probabilityInMEC) {
        if (numChainedMECs <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The number of added MECs per initial state must be at least 1");
        }
        else if (chainLengthInMEC <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The length of added MEC-Chains per initial state must be at least 1");
        }
        else if (sinkProbabilityPerMEC < 0 || sinkProbabilityPerMEC > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability of reaching a sink was set to "+sinkProbabilityPerMEC+" which is not in [0,1]");
        }
        else if (probabilityInMEC < 0 || probabilityInMEC > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability to branch in chain was set to "+probabilityInMEC+" which is not in [0,1]");
        }

        int numChains = 2;
        int numNewStatesPerMEC = 2 * numChainedMECs+1;
        int numNewStates = numChainedMECs*numNewStatesPerMEC;
        int addedSTPGInitialState = stpg.getNumStates();
        int currentInitialState = addedSTPGInitialState;
        int nextInitialState;

        int currentState;

        stpg.addStates(numNewStates);
        if (remain != null) remain.set(stpg.getNumStates()-numNewStates, stpg.getNumStates());
        addToStateList(stpg, stpg.getNumStates()-numNewStates, stpg.getNumStates());

        for (int chainedMEC = 0; chainedMEC < numChainedMECs; chainedMEC++) {
            // Make sure to know which will be the next initialState
            nextInitialState = currentInitialState + numNewStatesPerMEC;

            // Index out of bounds -> reached the end of prepended MEC. Next state the MEC leads into must be old initial state
            if (nextInitialState == stpg.getNumStates()) {
                nextInitialState = initialState;
            }

            stpg.setPlayer(currentInitialState, 2); //Min-Player gets to decide the beginning of MEC

            // Add transitions from chains initial state to both chains
            Distribution choiceDistribution = new Distribution();
            choiceDistribution.add(currentInitialState+1, 1.0);
            stpg.addChoice(currentInitialState, choiceDistribution);

            choiceDistribution = new Distribution();
            choiceDistribution.add(currentInitialState+2, 1.0);
            stpg.addChoice(currentInitialState, choiceDistribution);

            // Loop for the chain-complexes
            for (int chainMember = 0; chainMember < chainLengthInMEC; chainMember++) {
                // Loop for upper and lower chain
                for (int chain = 1; chain < numChains + 1; chain++) {
                    currentState = currentInitialState + chainMember * numChains + chain;


                    stpg.setPlayer(currentState, 1);
                    int previousState = currentState - numChains;
                    int nextState = currentState + numChains;

                    if (previousState <= currentInitialState) {
                        previousState = currentInitialState;
                    }
                    if (nextState >= currentInitialState + numNewStatesPerMEC) {
                        nextState = nextInitialState;
                    }

                    // add Choice to go back
                    choiceDistribution = new Distribution();
                    choiceDistribution.add(previousState, 1);
                    stpg.addChoice(currentState, choiceDistribution);

                    // add Choice to go to next
                    if (nextState == nextInitialState) {
                        // exit choices of the MEC
                        double probabilityNotToSink = 1.0 - sinkProbabilityPerMEC;
                        double probabilityToReachNextMEC = probabilityNotToSink * 0.5;
                        double probabilityToGoBackToCurrentInitialState = probabilityNotToSink * 0.5;
                        if (chain == 1) {
                            probabilityToReachNextMEC = probabilityNotToSink * 0.4;
                            probabilityToGoBackToCurrentInitialState = probabilityNotToSink * 0.6;
                        }
                        choiceDistribution = new Distribution();
                        if (sinkProbabilityPerMEC > 0) choiceDistribution.add(sinkState, sinkProbabilityPerMEC);
                        if (sinkProbabilityPerMEC < 1) {
                            choiceDistribution.add(nextInitialState, probabilityToReachNextMEC);
                            choiceDistribution.add(currentInitialState, probabilityToGoBackToCurrentInitialState);
                        }
                    }
                    else {
                        // If not exit, one can move forward or jump to other chain
                        int parallelState = (chain == 1) ? (currentState + 1) : (currentState - 1);
                        choiceDistribution = new Distribution();
                        choiceDistribution.add(nextState, 1.0 - probabilityInMEC);
                        if (probabilityInMEC > 0) choiceDistribution.add(parallelState, probabilityInMEC);
                    }
                    stpg.addChoice(currentState, choiceDistribution);
                }
            }
            currentInitialState = nextInitialState;
        }

        ModelExtensionResult extensionResult = new ModelExtensionResult(addedSTPGInitialState, stpg);
        return extensionResult;
    }

    private static ModelExtensionResult addStartingProbabilityGraph(STPGExplicit stpg, BitSet remain, int oldInitialState, int numComponents, int treeDepth, int treeBranchingFactor,
                                                                    double probabilityToReachComponentsInitialState, double probabilityToReachSink, int sinkState) {
        if (numComponents <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The number of added Components per initial state must be at least 1");
        }
        else if (treeDepth < 0) {
            throw new IllegalArgumentException("[Model Extension]: The depth of added ProbabilityTrees per initial state must be at least 0");
        }
        else if (treeBranchingFactor <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The treeBranchingFactor must be at least 1");
        }
        else if (probabilityToReachSink < 0 || probabilityToReachSink > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability of reaching a sink was set to "+probabilityToReachSink+" which is not in [0,1]");
        }
        else if (probabilityToReachComponentsInitialState < 0 || probabilityToReachComponentsInitialState > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability to branch in tree was set to "+probabilityToReachComponentsInitialState+" which is not in [0,1]");
        }

        Random random = new Random();
        random.setSeed(100); //Set deterministic Seed because we do want deterministic models but don't care about who the states belong

        int numNewStatesPerComponent = 1;
        for (int depth = 1; depth <= treeDepth; depth++) {
            numNewStatesPerComponent += Math.pow(treeBranchingFactor, depth);
        }
        int numNewStates = numNewStatesPerComponent*numComponents;
        int addedSTPGInitialState = stpg.getNumStates();
        int numStatesBeforeExtension = stpg.getNumStates();
        int currentInitialState = addedSTPGInitialState;
        int nextInitialState;

        int currentState;

        stpg.addStates(numNewStates);
        if (remain != null) remain.set(stpg.getNumStates()-numNewStates, stpg.getNumStates());
        addToStateList(stpg, stpg.getNumStates()-numNewStates, stpg.getNumStates());

        for (int component = 0; component < numComponents; component++) {
            // Make sure to know which will be the next initialState
            nextInitialState = currentInitialState + numNewStatesPerComponent;

            // Index out of bounds -> reached the end of prepended MEC. Next state the MEC leads into must be old initial state
            if (nextInitialState == stpg.getNumStates()) {
                nextInitialState = oldInitialState;
            }

            for (currentState = currentInitialState; currentState < currentInitialState + numNewStatesPerComponent; currentState++) {
                stpg.setPlayer(currentState, random.nextInt(2)+1); //Set player to either 1 or 2
                int childrenOffset = (currentState - currentInitialState) * treeBranchingFactor + currentInitialState; //Indexing in Trees
                for (int nextState = childrenOffset + 1; nextState <= childrenOffset + treeBranchingFactor; nextState++) {
                    if (nextState >= currentInitialState + numNewStatesPerComponent) {
                        nextState = nextInitialState;
                    }
                    Distribution transitionDistribution = new Distribution();

                    double probabilityNotToSink = 1.0 - probabilityToReachSink;
                    if (probabilityToReachSink > 0) {
                        transitionDistribution.add(sinkState, probabilityToReachSink);
                    }
                    if (probabilityToReachSink < 1) {
                        transitionDistribution.add(currentInitialState, probabilityNotToSink * probabilityToReachComponentsInitialState);
                        transitionDistribution.add(nextState, probabilityNotToSink * (1-probabilityToReachComponentsInitialState));
                    }
                    stpg.addChoice(currentState, transitionDistribution);

                    //Without this statement the for-loop may never end
                    if (nextState == nextInitialState) {
                        break;
                    }
                }
            }
            currentInitialState = nextInitialState;
        }

        ModelExtensionResult extensionResult = new ModelExtensionResult(addedSTPGInitialState, stpg);
        return extensionResult;
    }

    private static void addToStateList(STPGExplicit stpg, int statesFrom, int statesTo) {
        for (int i = statesFrom; i<=statesTo; i++) {
            State s = new State(1);
            s.setValue(0, (Integer) i);
            stpg.statesList.add(s);
        }
    }

    public static STPG createCopyOfExistingSTPG(STPG stpg) {
        STPGExplicit newSTPG = new STPGExplicit();

        newSTPG.addStates(stpg.getNumStates());

        for (int state = 0; state < stpg.getNumStates(); state++) {
            newSTPG.setPlayer(state, stpg.getPlayer(state));
            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                Distribution distribution = getDistribution(stpg, state, choice);
                newSTPG.addChoice(state, distribution);
            }
        }

        for (int initialState : stpg.getInitialStates()) {
            newSTPG.addInitialState(initialState);
        }
        return newSTPG;
    }



    public static Distribution getDistribution(STPG stpg, int state, int action) {
        Distribution d = new Distribution();
        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, action); it.hasNext(); ) {
            Map.Entry<Integer, Double> tr = it.next();
            d.add(tr.getKey(), tr.getValue());
        }
        return d;
    }

    private static class ModelExtensionResult {
        public int addedInitialState;
        public STPGExplicit resultingSTPG;

        public ModelExtensionResult(int addedInitialState, STPGExplicit resultingSTPG) {
            this.addedInitialState = addedInitialState;
            this.resultingSTPG = resultingSTPG;
        }
    }

    private static class ModelExtensionParameters {
        public boolean extendWithMECs;
        public boolean extendWithProbabilityComponents;

        //MEC Parameters
        int mecExtension_numMECs;
        int mecExtension_chainLength;
        double mecExtension_probabilityLeadingToSink;
        double mecExtension_probabilityToGoBackToIntialState;
        double mecExtension_probabilityOfChainSwitch;

        //ProbabilityTree Parameters
        int probExtension_numComponents;
        int probExtension_componentBranchingFactor;
        int probExtension_componentTreeDepth;
        double probExtension_probabilityLeadingToSink;
        double probExtension_probabilityToGoBackToIntialState;
    }

    private static class ConfigParserJSON
    {
        @SuppressWarnings("unchecked")
        public static ModelExtensionParameters read(String filePath) throws PrismException
        {
            ModelExtensionParameters config = new ModelExtensionParameters();
            //JSON parser object to parse read file
            JSONParser jsonParser = new JSONParser();

            File configFile = null;
            try {
                configFile = new File(filePath);
                if (!configFile.exists()) {
                    throw new PrismException("[Model Extension]: Provided filepath "+filePath+" is not a File");
                }
            }
            catch (SecurityException e) {
                e.printStackTrace();
            }

            try (FileReader reader = new FileReader(filePath))
            {
                //Read JSON file
                Object obj = jsonParser.parse(reader);

                JSONObject extensions = (JSONObject) obj;

                JSONObject extensionMEC = (JSONObject) extensions.get("MECModel");

                config.extendWithMECs = Boolean.parseBoolean(extensionMEC.get("use").toString());
                config.mecExtension_numMECs = Integer.parseInt(extensionMEC.get("numMECs").toString());
                config.mecExtension_chainLength = Integer.parseInt(extensionMEC.get("chainLength").toString());
                config.mecExtension_probabilityLeadingToSink = Double.parseDouble(extensionMEC.get("probabilityLeadingToSink").toString());
                config.mecExtension_probabilityToGoBackToIntialState = Double.parseDouble(extensionMEC.get("probabilityToGoBackToInitialState").toString());
                config.mecExtension_probabilityOfChainSwitch = Double.parseDouble(extensionMEC.get("probabilityForChainSwitch").toString());

                JSONObject extensionProbabilistic = (JSONObject) extensions.get("ProbabilisticModel");

                config.extendWithProbabilityComponents = Boolean.parseBoolean(extensionProbabilistic.get("use").toString());
                config.probExtension_numComponents = Integer.parseInt(extensionProbabilistic.get("numComponents").toString());
                config.probExtension_componentBranchingFactor = Integer.parseInt(extensionProbabilistic.get("componentBranchingFactor").toString());
                config.probExtension_componentTreeDepth = Integer.parseInt(extensionProbabilistic.get("componentTreeDepth").toString());
                config.probExtension_probabilityLeadingToSink = Double.parseDouble(extensionProbabilistic.get("probabilityLeadingToSink").toString());
                config.probExtension_probabilityToGoBackToIntialState = Double.parseDouble(extensionProbabilistic.get("probabilityToGoBackToInitialState").toString());

            } catch (Exception e) {
                e.printStackTrace();
                throw new PrismException("Reading the Model Extension File failed");
            }
            return config;
        }
    }

    //Taken from: https://mkyong.com/java/how-to-read-xml-file-in-java-dom-parser/
    private static class ConfigParserXML {


        public static ModelExtensionParameters read(String filepath) throws PrismException {
            ModelExtensionParameters config = new ModelExtensionParameters();

            // Instantiate the Factory
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            File configFile = null;
            try {
                configFile = new File(filepath);
                if (!configFile.exists()) {
                    throw new PrismException("[Model Extension]: Provided filepath "+filepath+" is not a File");
                }
            }
            catch (SecurityException e) {
                e.printStackTrace();
            }

            try {

                // optional, but recommended
                // process XML securely, avoid attacks like XML External Entities (XXE)
                dbf.setFeature(XMLConstants.FEATURE_SECURE_PROCESSING, true);

                // parse XML file
                DocumentBuilder db = dbf.newDocumentBuilder();

                Document doc = db.parse(configFile);

                // optional, but recommended
                // http://stackoverflow.com/questions/13786607/normalization-in-dom-parsing-with-java-how-does-it-work
                doc.getDocumentElement().normalize();

                // Get MEC config
                NodeList list = doc.getElementsByTagName("MECModel");
                for (int temp = 0; temp < list.getLength(); temp++) {
                    Node node = list.item(temp);
                    if (node.getNodeType() == Node.ELEMENT_NODE) {
                        Element element = (Element) node;

                        config.extendWithMECs = Boolean.parseBoolean(element.getElementsByTagName("use").item(0).getTextContent());
                        config.mecExtension_numMECs = Integer.parseInt(element.getElementsByTagName("numMECs").item(0).getTextContent());
                        config.mecExtension_chainLength = Integer.parseInt(element.getElementsByTagName("chainLength").item(0).getTextContent());
                        config.mecExtension_probabilityLeadingToSink = Double.parseDouble(element.getElementsByTagName("probabilityLeadingToSink").item(0).getTextContent());
                        config.mecExtension_probabilityToGoBackToIntialState = Double.parseDouble(element.getElementsByTagName("probabilityToGoBackToInitialState").item(0).getTextContent());
                        config.mecExtension_probabilityOfChainSwitch = Double.parseDouble(element.getElementsByTagName("probabilityForChainSwitch").item(0).getTextContent());

                    }
                }

                // Get Probabilistic config
                list = doc.getElementsByTagName("ProbabilisticModel");
                for (int temp = 0; temp < list.getLength(); temp++) {
                    Node node = list.item(temp);
                    if (node.getNodeType() == Node.ELEMENT_NODE) {
                        Element element = (Element) node;

                        config.extendWithProbabilityComponents = Boolean.parseBoolean(element.getElementsByTagName("use").item(0).getTextContent());
                        config.probExtension_numComponents = Integer.parseInt(element.getElementsByTagName("numComponents").item(0).getTextContent());
                        config.probExtension_componentBranchingFactor = Integer.parseInt(element.getElementsByTagName("componentBranchingFactor").item(0).getTextContent());
                        config.probExtension_componentTreeDepth = Integer.parseInt(element.getElementsByTagName("componentTreeDepth").item(0).getTextContent());
                        config.probExtension_probabilityLeadingToSink = Double.parseDouble(element.getElementsByTagName("probabilityLeadingToSink").item(0).getTextContent());
                        config.probExtension_probabilityToGoBackToIntialState = Double.parseDouble(element.getElementsByTagName("probabilityToGoBackToInitialState").item(0).getTextContent());

                    }
                }

            } catch (ParserConfigurationException | SAXException | IOException e) {
                e.printStackTrace();
            }
            return config;
        }
    }
}