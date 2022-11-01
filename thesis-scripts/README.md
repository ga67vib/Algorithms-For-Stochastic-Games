

Usually, to each result folder there exists a python file with the same name.
The python file contains the timeout and memory settings and provides information about which models where used:
- `random_input_dir` contains the path to the random models that were used
    - Random models are only used if the random model code block in around line 165 or 175 is not commented out
- `extension_config_folder_path` contains the path to the folder of configs used to create the models that are extended during Java runtime
    - Configs are only used if the code block around the `extension_config_folder_path` is not commented out
- Paths to handcrafted and real models are hardcoded in the python file. Thus they are used iff they are not commented out

### Result Folders

Next, we provide a short overview of the result folders. Result folders that we have used in publications are marked in bold.
[A] and [M] indicate in which publication we have used which folder.

[A]: Atva Paper - "Optimistic and Topological Value Iteration for Simple Stochastic Games"
[M]: Alexander Slivinskiys Master Thesis - "Facilitating Experimental Analysis of Algorithms for Stochastic Games: Random Model Generation and Mining the Results"

- **bigConfModels**: simple tree-like models with many states that are grouped either into one or multiple SCCs
    - [A] Used in Table 1, Figs. 3b, 
    - [M] Used in Figures 7.6a, 7.6b, 7.7, 7.8, 7.9, 7.13b
- ~bigRealModels: real models with parameter settings such that the models are considered large and hard to compute
- **handcraftedModels**: handcrafted models
    - [A] Used in Figs. 3a, 4a, 7, 8, 9, 10a
    - [M] Used in Figs 7.5, 7.6a, 7.7, 7.8, 7.9, 7.10, 7.11, 7.12, 7.13a, 7.14
- ~midConfModels2: like bigConfModels, but with smaller state spaces
- **oviGood**: Handcrafted edge case where OVI is better than BVI (but not better than WP3)
    - [A] Used in Fig. 11
- **randomBigSCC**: Randomly generated models with one big SCC
    - [M] Used in Fig. 7.7, 7.8, 7.9
- **randomRandomThesis**: Simple randomly generated models. Are solved fast by everyone usually
    - [A] Used in Figs 7, 8, 9, 10b
    - [M] Used in Figs 7.2, 7.4b, 7.5, 7.6a, 7.7, 7.8, 7.9, 7.10, 7.11, 7.12, 7.13a, 7.14
- **randomRandom10Act**: Randomly generated models where every state has around 10 actions. More complicated to solve than randomRandomThesis
    - [A] Used in Figs. 2b, 3a, 4b, 6, 7, 8, 9, 10b
    - [M] Used in Figs 7.4b, 7.5, 7.6a, 7.7, 7.8, 7.9, 7.10, 7.11, 7.12, 7.13a, 7.14
- **randomSCCFixed2**: Randomly generated models with the scc guideline
    - [M] Used in Figs 7.3, 7.4b, 7.5, 7.6a, 7.7, 7.8, 7.9, 7.10, 7.11, 7.12, 7.13a, 7.14
- randomTreeThesis: Randomly generated models with the tree guideline
- **realModels**: Real models
    - [A] Used in Figs. 2a, 3a, 4a, 5, 7, 8, 9, 10a
    - [M] Used in Figs. 7.1, 7.4a, 7.5, 7.6a, 7.7, 7.8, 7.9, 7.10, 7.11, 7.12, 7.13a, 7.14