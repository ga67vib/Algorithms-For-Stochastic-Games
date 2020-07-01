# Algorithms-For-Stochastic-Games

Extension of PRISM-games (link).
Algorithms as described in GandALF submission "Title".

## Installation

Standard PRISM tutorial.
For QP stuff [TODO: Sasha write tutorial]
<br/>[Not sure where to place that we expect the people (also for this guide) to have Linux]
### Using quadratic programming
To be able to use our quadratic programming implementation, you need a Gurobi 9.0.0 or more recent license file. We use Gurobi to construct the quadratic program from the stochastic game as well as to solve it.
To obtain a license, check out the following website:<br/>
https://www.gurobi.com/documentation/9.0/quickstart_linux/obtaining_a_grb_license.html#section:ObtainLicense<br/>
Gurobi provides free academic licenses that expire after 1 year.

Follow the guidelines on https://www.gurobi.com/documentation/quickstart.html to get Gurobi, to obtain a license and to activate it using `grbgetkey`.
Now you should be able to execute following command from the folder prism-games-3.0.beta-src/prism to check whether everything works:<br/>
`./bin/prism ../../case_studies/BigMec.prism ../../case_studies/BigMec.props -const N=1 -qp`

[TODO: LD_LIBRARY_PATH, INSTALL GUROBI?]

CPLEX:
[ACTUALLY I'M NOT SURE HOW TO SHIP CPLEX, SINCE THEY DON'T REALLY REQUIRE A LICENSE]

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
