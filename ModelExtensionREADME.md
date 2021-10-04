#Using Model Extensions

To use a model extension on your SMG, add the switch -smg_extend <PATH_TO_CONFIG> to the execution command and provide the relative path to the config.json that contains the extension parameters

#Adding new Extensions

To add a new Extension Model, you have to do the following steps:
1. Add a new class "SMGModelExtension_<Your Model Name>" to src/explicit/smgModelExtensions that extends SMGModelExtension
2. The class should get all its necessary parameters in the constructor, and must do the model extension in the method "extendSMG"
	- Don't forget to extend the statesList and the remain Bitset
	- You should not override the initial states of the STPG in "extendSMG"
3. Next, we extend the JSON-config to include your new parameters:
	- Extend the class "ModelExtensionParameters" by the parameters you need for your Model
	- Extend the config file. Do not delete the other model parameters. If you do not want to use the other extensions, simply set "use" of these models to "false"
	- Extend the method read of the ConfigParserJSON class in the SMGModelExtender file
4. Add a new instance of your Model to the array "extensions" in the method "extendModel" of class "SMGModelExtender"
