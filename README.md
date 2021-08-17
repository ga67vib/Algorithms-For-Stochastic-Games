# Algorithms-For-Stochastic-Games

Extension of the project PRISM-games [https://github.com/prismmodelchecker/prism-games] with the algorithms as described in GandALF'20 paper "Comparison of Algorithms for Simple Stochastic Games" and the FSTTCS'21 submission "Faster Value Iteration for Simple Stochastic Games".

## License:

For our modifications of the code, we give an MIT license. For the code that was originally part of PRISM-games, we refer to their license file in prism-games-3.0.beta-src/COPYING.txt and the "Licensing" section of their README in prism-games-3.0.beta-src/README.md.

## Installation

### Dependencies
- Java 8 or more recent
- Python 3.7 or more recent
- C and C++ compiler
- PRISM-Games (https://www.prismmodelchecker.org/games/)
  - PPL
- Gurobi Optimizer 9.0.0 or more recent (Optional) (https://www.gurobi.com/)
- CPLEX (Optional) (https://www.ibm.com/products/ilog-cplex-optimization-studio)
- AMPL (Optional) (https://ampl.com/)

**We assume that you are using a Linux distribution.** For Windows and Mac OS users the guides we refer to also provide instructions, but you may encounter difficulties putting everything together.

### Setting up PRISM-games

Our implementation is an extension of PRISM-games.
Thus, we require PPL and PRISM to be installed.
Follow the guidelines given at:<br/>
https://www.prismmodelchecker.org/games/installation.php<br/>
Now the following command should execute without errors from the prism-games-3.0.beta-src/prism folder:
`./bin/prism`
At this point you should be able to run the case-study script for everything except for the mathematical programming approaches.

### More tutorial on setting up PRISM-games
- Download PRISM sources (not binaries), e.g. from https://github.com/ga67vib/Algorithms-For-Stochastic-Games or https://prismmodelchecker.org/games/

- Install PRISM (more instructions available here: https://prismmodelchecker.org/games/installation.php)
    - Open project in IntelliJ (use prism folder inside prism-games-3-sthsth-whatever folder). Note: This is important. If you select choose right, IntelliJ finds some libraries and stuff itself. Otherwise you'll have to manually add them or tell it where sources are and stuff.
    - From prism folder, run make
        Possibly encounter weird errors. Fix by googling or asking supervisor
    - Check that installation worked: execute (from prism folder) bin/prism in terminal. Should report version and stuff. Congrats, you can now run prism!

- How to be able to run it from IntelliJ:
    - Build the project (There's a menu at the top)
        Possible errors we know:<br>
            a: Rightclick on the src-folder -> Mark Directory as -> Sources root<br>
            b: Project settings (Ctrl+Alt+Shift+S) -> Set output directory as prism/classes<br>
            c: Also in project settings, select correct JDK (and make sure to use java 8, at least if using prism-games)<br>
            d: Also in project settings, add lib folder to libraries<br>
    - Add a configuration to run it, i.e. top right there should be "Add configuration" or "edit configuration". Add new configuration of type "Application". Set the Main class to prism.PrismCL and add a new Environment variable - LD_LIBRARY_PATH=:lib. Note that on MAC OS, this is called DYLD_LIBRARY_PATH.<br>
            In case of weird errors, check a, b, d from before (often stuff with libraries missing, so copy them from some other place. Or source folder not specified (can't select PrismCL as main class), so select it. Or wrong Java/no SDK selected)<br>
            In case of errors that include "Parma Polyhedra Library", you still have to properly install the ppl library, which is needed for multi-objective stuff (and hence sadly also in general, cause we have to compile all of PRISM)<br>
                In the Algorithms-for-SG github, PPL is already in ext folder, so it might also work out of the box (and break if you try to install it manually)<br>
            If there are libraries missing (PPL or commons or sth like that), then go to project settings (Ctrl Alt Shift S) -> Libraries -> + and add the ext and lib folder. Or PPl.jar directly. 
    - Click play and see if it runs. If yes: Party. If no: Call someone who might know what to do, e.g. Pranav, Tobi or Maxi.

- Now to add command line parameters in order to run the specific thing you want
    Syntax when calling from command line: bin/prism <path/to/model> <path/to/property> -const <constants> <configuration modifiers>
    In IntelliJ, you add everything but the "bin/prism" in front to the "command line arguments" field of the configuration.
    E.g. put in 
        ../../case-studies/BigMec.prism ../../case-studies/BigMec.props -const N=1 -ii -smg_opts 2
    in order run on the model BigMec with the BigMec.props property file. We set the free scaling constant of the model N to 1. We use -ii, so bounded value iteration (aka interval iteration). And we set smg_opts to 2, which the method can use to differentiate variants of the algorithm. E.g. 2 corresponds to sound value iteration (unless we changed it in the meantime).

- How to modify the code
    Find the method that actually does the stuff you are interested in. Often this is in explicit/STPGModelChecker (at least for my research), as this class contains the algorithms to model check an STPG. 
    
- How to pass arguments through PRISM
    Breakpoint in PrismCL.main, debug, pray


### Using quadratic programming
To be able to use our quadratic programming implementation, you need a Gurobi 9.0.0 or more recent license file. We use Gurobi to construct the quadratic program from the stochastic game as well as to solve it. Gurobi provides free academic licenses that expire after one year.
To obtain a license, install gurobi and activate the license follow the guidelines provided on:<br/>
https://www.gurobi.com/documentation/quickstart.html<br/>
After executing `./grbgetkey <YOUR_LICENSE_KEY>` in the bin folder of Gurobi, you should be able to use quadratic programming to solve simple stochastic games with our implementation. To verify this, try out the following command in the prism-games-3.0.beta-src/prism folder:<br/>
`./bin/prism ../../case_studies/BigMec.prism ../../case_studies/BigMec.props -const N=1 -qp`

#### CPLEX
Due to license requirements, we disabled solving simple stochastic games with CPLEX.
If you want to use CPLEX to solve simple stochastic games with quadratic programming, follow these steps:
- Download IBM ILOG CPLEX Optimization Studio (They provide free academic licences): https://www.ibm.com/products/ilog-cplex-optimization-studio
- Use `sudo chmod +x cplex_studioZZZ.YYY.bin` to make the installer executable (ZZZ represents the version, YYY your system)
- Execute the .bin to install CPLEX
- Copy the following two files into the prism-games-3.0.beta-src/prism/lib folder:
  - cplex.jar which can be found in cplex/lib
  - libcplex12100.so which can be found in cpelx/bin/<YourSystem>
- Remove the comments that are marked with "CPLEX-RELATED" in the file SMGPolyProgSolverGurobi.java in prism-games-3.0.beta-src/prism/src/explicit
  - Exchange the `throwCplexError()`-calls in the file by the method `solveQPCplex()`
- If you have already built PRISM, you will have to use `make clean && make` in prism-games-3.0.beta-src/prism
Once you have followed these steps, the following command should work without errors in the prism-games-3.0.beta-src/prism folder:<br/>
`./bin/prism ../../case_studies/BigMec.prism ../../case_studies/BigMec.props -const N=1 -qp -smg_opts 3`

### Using higher order programming
To be able to use higher order programming, you require an AMPL license. We use AMPL to formulate the stochastic game as a nonlinear program and Minos to compute the optimal solution.
To obtain a license, check out the following website:<br/>
https://ampl.com/products/ampl/ampl-for-research/<br/>
AMPL provides a free 30-day trail-license and a student license. Once you have obtained the license, put it into prism-games-3.0.beta-src/prism. In this folder you should be able to get information about your license with:<br/>
`./ampl -v`<br/>
Note that we expect you to use prism for AMPL from either the prism-games-3.0.beta-src/prism folder (via `./bin/prism`) or from the scripts-folder using for example the run_benchmarks.py script. Otherwise you may encounter difficulties since the program won't be able to find its path to the AMPL executable.

## Running the code

You can either use the script scripts/run_benchmarks.py to perform multiple (or all) experiments at once, see scripts/README.md for more details.
Or you can directly execute PRISM from the command line. The necessary switches are -ii for BVI, -qp for mathematical programming and -politer for strategy iteration.
The switch -smg_opts can be used to select the optimizations. Refer to the runscript to see which number corresponds to which combination of optimizations.

### Pueue

Download pueue (follow instructions here https://github.com/Nukesor/pueue) for the fsttcs runscript. It is an awesome utility for parallelizing commands.

Then execute '''pueued -d''', which starts the pueue daemon and makes it listen for tasks.
Configure pueue as described here: https://github.com/Nukesor/pueue/wiki/Configuration and here https://github.com/Nukesor/pueue/wiki/Groups, setting the number of parallel tasks and using groups to ensure every task has a single CPU.

Then run the script, which adds the task to the pueue and then runs them in parellel.

## Looking at the code

Our new algorithms are located in prism-games-3.0.beta-src/prism/src/explicit/STPGModelChecker.java, namely computeReachProbsValIter (with bounded=true), computeReachProbsPolIter and computeReachProbsQuadProg.
