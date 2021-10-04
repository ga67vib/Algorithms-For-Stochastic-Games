package explicit;

import java.io.FileReader;
import java.util.*;

import explicit.Distribution;
import explicit.STPG;
import explicit.STPGExplicit;
import explicit.smgModelExtensions.*;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
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

        int sink = newStpg.addState();
        Distribution sinkSelfLoop = new Distribution();
        sinkSelfLoop.add(sink, 1.0);
        newStpg.addChoice(sink, sinkSelfLoop);
        if (remain != null) remain.set(sink);
        addToStateList(newStpg, sink, sink);

        Iterable<Integer> initialStates;

        List<State> statesList = ((STPGExplicit) stpg).statesList;
        SMGModelExtension[] extensions = {
                new SMGModelExtension_BigMEC(
                        newStpg,
                        remain,
                        statesList,
                        config.extendWithBigMECs,
                        config.mecExtension_numMECs,
                        config.mecExtension_chainLength,
                        config.mecExtension_probabilityLeadingToSink,
                        sink,
                        config.mecExtension_probabilityOfChainSwitch),
                new SMGModelExtension_ProbabilityTree(
                        newStpg,
                        remain,
                        statesList,
                        config.extendWithProbabilityTrees,
                        config.probExtension_numComponents,
                        config.probExtension_componentTreeDepth,
                        config.probExtension_componentBranchingFactor,
                        config.probExtension_probabilityToGoBackToIntialState,
                        config.probExtension_probabilityLeadingToSink,
                        sink),
                new SMGModelExtension_ActionTree(
                        newStpg,
                        remain,
                        statesList,
                        config.extendWithActionTrees,
                        config.actionTreeExtension_numComponents,
                        config.actionTreeExtension_componentNumStates,
                        config.actionTreeExtension_componentBranchingFactor,
                        config.actionTreeExtension_probabilityToGoBackToIntialState,
                        config.actionTreeExtension_probabilityLeadingToSink,
                        sink)
        };

        ArrayList<Integer> newInitialStates = new ArrayList<>();
        for (SMGModelExtension modelExtension : extensions) {
            if (!modelExtension.shouldUseExtension()) continue;

            initialStates = newStpg.getInitialStates();
            for (int initialState : initialStates) {
                ModelExtensionResult modelExtensionResult = modelExtension.extendSMG(initialState);
                newStpg = modelExtensionResult.resultingSTPG;
                newInitialStates.add(modelExtensionResult.addedInitialState);
            }

            newStpg.clearInitialStates();
            for (int newInitialState : newInitialStates) {
                newStpg.addInitialState(newInitialState);
            }
            newInitialStates.clear();
        }

        return newStpg;
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

                config.extendWithBigMECs = Boolean.parseBoolean(extensionMEC.get("use").toString());
                config.mecExtension_numMECs = Integer.parseInt(extensionMEC.get("numMECs").toString());
                config.mecExtension_chainLength = Integer.parseInt(extensionMEC.get("chainLength").toString());
                config.mecExtension_probabilityLeadingToSink = Double.parseDouble(extensionMEC.get("probabilityLeadingToSink").toString());
                config.mecExtension_probabilityToGoBackToIntialState = Double.parseDouble(extensionMEC.get("probabilityToGoBackToInitialState").toString());
                config.mecExtension_probabilityOfChainSwitch = Double.parseDouble(extensionMEC.get("probabilityForChainSwitch").toString());

                JSONObject extensionProbabilistic = (JSONObject) extensions.get("ProbabilisticModel");

                config.extendWithProbabilityTrees = Boolean.parseBoolean(extensionProbabilistic.get("use").toString());
                config.probExtension_numComponents = Integer.parseInt(extensionProbabilistic.get("numComponents").toString());
                config.probExtension_componentBranchingFactor = Integer.parseInt(extensionProbabilistic.get("componentBranchingFactor").toString());
                config.probExtension_componentTreeDepth = Integer.parseInt(extensionProbabilistic.get("componentTreeDepth").toString());
                config.probExtension_probabilityLeadingToSink = Double.parseDouble(extensionProbabilistic.get("probabilityLeadingToSink").toString());
                config.probExtension_probabilityToGoBackToIntialState = Double.parseDouble(extensionProbabilistic.get("probabilityToGoBackToInitialState").toString());

                JSONObject extensionActionTree = (JSONObject) extensions.get("ActionTreeModel");

                config.extendWithActionTrees = Boolean.parseBoolean(extensionActionTree.get("use").toString());
                config.actionTreeExtension_numComponents = Integer.parseInt(extensionActionTree.get("numComponents").toString());
                config.actionTreeExtension_componentBranchingFactor = Integer.parseInt(extensionActionTree.get("componentBranchingFactor").toString());
                config.actionTreeExtension_componentNumStates = Integer.parseInt(extensionActionTree.get("componentNumberOfStates").toString());
                config.actionTreeExtension_probabilityLeadingToSink = Double.parseDouble(extensionActionTree.get("probabilityLeadingToSink").toString());
                config.actionTreeExtension_probabilityToGoBackToIntialState = Double.parseDouble(extensionActionTree.get("probabilityToGoBackToInitialState").toString());


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

                        config.extendWithBigMECs = Boolean.parseBoolean(element.getElementsByTagName("use").item(0).getTextContent());
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

                        config.extendWithProbabilityTrees = Boolean.parseBoolean(element.getElementsByTagName("use").item(0).getTextContent());
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