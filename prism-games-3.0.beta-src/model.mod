#variable declarations
var s0 >=0, <=1;
param s1 =1; #Yes
param s2 =0; #No


#Constraints for non-MEC, non-trivial states
subject to const_s0_a0: s0 >= 0.5 + 0;
subject to const_s0_a1: s0 <= 0.5 + 0;

#MEC-Constraints

#Objective function
minimize value: 0;
