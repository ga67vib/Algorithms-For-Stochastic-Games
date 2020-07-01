# Algorithms-For-Stochastic-Games

Extension of PRISM-games (link).
Algorithms as described in GandALF submission "Title".

## Installation

### Dependencies
- Java 8 or more recent
- Python 3.7 or more recent
- C and C++ compiler
- PRISM-Games (https://prismmodelchecker.org/games/installation.php)
  - PIL
- Gurobi Optimizer 9.0.0 or more recent (Optional) (https://www.gurobi.com/)
- AMPL (Optional) (https://ampl.com/)

**We assume that you are using a Linux distribution.** For Windows and Mac OS users the guides we refer to also provide instructions, but you may encounter difficulties putting everything together.

### Setting up PRISM-games

Our implementation is an extension of PRISM-games.
Thus, we require PPL to be installed.
Follow the guidelines given at:<br/>
https://prismmodelchecker.org/games/installation.php<br/>
Now the following command should execute without errors from the prism-games-3.0.beta-src/prism folder:
`./bin/prism`
At this point you should be able to run the case-study script for everything except for the mathematical programming approaches.


### Using quadratic programming
To be able to use our quadratic programming implementation, you need a Gurobi 9.0.0 or more recent license file. We use Gurobi to construct the quadratic program from the stochastic game as well as to solve it. Gurobi provides free academic licenses that expire after one year.
To obtain a license, install gurobi and activate the license follow the guidelines provided on:<br/>
https://www.gurobi.com/documentation/quickstart.html<br/>
After executing `./grbgetkey <YOUR_LICENSE_KEY>` in the bin folder of Gurobi, you should be able to use quadratic programming to solve simple stochastic games with our implementation. To verify this, try out the following command in the prism-games-3.0.beta-src/prism folder:<br/>
`./bin/prism ../../case_studies/BigMec.prism ../../case_studies/BigMec.props -const N=1 -qp`

####CPLEX
Due to license requirements, we disabled solving simple stochastic games with CPLEX

### Using higher order programming
To be able to use higher order programming, you will require an AMPL license. We use AMPL to formulate the stochastic game as a nonlinear program and Minos to compute the optimal solution.
To obtain a license, check out the following website:<br/>
https://ampl.com/products/ampl/ampl-for-research/<br/>
AMPL provides a free 30-day trail-license and a student license. Once you have obtained the license, put it into prism-games-3.0.beta-src/prism. In this folder you should be able to get information about your license with:<br/>
`./ampl version`<br/>
Note that we expect you to use prism for AMPL from either the prism-games-3.0.beta-src/prism folder (via `./bin/prism`) or from the scripts-folder using for example the run_benchmarks.py script. Otherwise you may encounter difficulties since the program won't be able to find it's path to the AMPL executable.

## TODO:
- put code here (Merging master and finish by hand as appropriate, excluding class files and the like)
- Check license of PRISM, check where to add our names
- Fix things in Readme
- Try it out
- If time permits: Make it work even without license; maybe by putting a dummy license file called 'your-license-here' and telling gurobi/AMPL to look there
- scripts and case studies
- README should tell people where they find our stuff
- polish code

## License:

For our modifications of the code, we have give an MIT license. For the code that was originally part of PRISM-games, we refer to their license file in prism-games-3.0.beta-src/COPYING.txt and the "Licensing" section of their README in prism-games-3.0.beta-src/COPYING.txtprism-games-3.0.beta-src/README.md
