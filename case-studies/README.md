These case studies are either from PRISM-games (http://www.prismmodelchecker.org/games/casestudies.php), from the paper [https://dblp.uni-trier.de/rec/bibtex/journals/corr/abs-2005-04018] or handcrafted for the GandALF paper.

More concretely:

The case studies coins, prison_dil, adt, charlton, cdmsn, cloud, mdsm, dice and two_investors are distributed with PRISM-games 3 or available on their case-study-website [http://www.prismmodelchecker.org/games/casestudies.php].

To judge the impact of a single small MEC, we prepended dice and cdmsn with a single MEC. The exits lead to the initial state of the original model with some probability, and the remaining probability leads to a sink. 

HW and AV are the models used in [https://dblp.uni-trier.de/rec/bibtex/journals/corr/abs-2005-04018]; the first two indices show the size of the grid, the last index denotes the single property used.

As interesting handcrafted examples, we used the adversarial model for value iteration from [https://dblp.uni-trier.de/rec/bibtex/journals/tcs/HaddadM18] (called hm) as well as two newly handcrafted models with either one large MEC (BigMec) or many 3 state MECs (MulMECs).
In BigMec, there is a MEC with two chains of N Maximizer states.
In MulMec, a single MEC is repeated N times. 
