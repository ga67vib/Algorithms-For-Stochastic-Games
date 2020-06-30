# Algorithms-For-Stochastic-Games

Extension of PRISM-games (link).
Algorithms as described in GandALF submission "Title".

## Installation

Standard PRISM tutorial.
For QP stuff [TODO: Sasha write tutorial]
### Using quadratic programming
To be able to use our quadratic programming implementation, you need a Gurobi 9.0.0 or more recent license file. We use Gurobi to construct the quadratic program from the stochastic game as well as to solve it.
To obtain a license, check out the following website:
https://www.gurobi.com/documentation/9.0/quickstart_linux/obtaining_a_grb_license.html#section:ObtainLicense
Gurobi provides free academic licenses that expire after 1 year.

Follow the guidelines on https://www.gurobi.com/documentation/quickstart.html to get Gurobi, to obtain a license and to activate it using grbgetkey.
Now you should be able to execute following command to check whether everything works:
[SOME EXAMPLE TEST-LINE like ./bin/prism ../prism-examples/smgs/simple/dice.prism ../prism-examples/smgs/simple/dice.props -const N=1 -qp]

[TODO: LD_LIBRARY_PATH]

CPLEX:
[ACTUALLY I'M NOT SURE HOW TO SHIP CPLEX, SINCE THEY DON'T REALLY REQUIRE A LICENSE]

### Using higher order programming
To be able to use higher order programming, you will require a AMPL license. We use AMPL to formulate the stochastic game as a nonlinear program and Minos to compute the optimal solution.
To obtain a license, check out the following website:
https://ampl.com/products/ampl/ampl-for-research/
AMPL provides a free 30-day trail-license and a student license.
[Check if path-stuff is all right]

## TODO:
- put code here (Merging master and finish by hand as appropriate, excluding class files and the like)
- Check license of PRISM, check where to add our names
- Fix things in Readme
- Try it out
- If time permits: Make it work even without license; maybe by putting a dummy license file called 'your-license-here' and telling gurobi/AMPL to look there
- scripts and case studies
- README should tell people where they find our stuff
- polish code
