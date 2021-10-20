#variable declarations
var s0 >=0, <=1;
param s1 =1; #Yes
var s2 >=0, <=1;
var s3 >=0, <=1;
param s4 =1; #Yes
param s5 =0; #No
param s6 =0; #No
param s7 =1; #Yes
var s8 >=0, <=1;
param s9 =1; #Yes
param s10 =1; #Yes
param s11 =0; #No
param s12 =1; #Yes
param s13 =0; #No
param s14 =1; #Yes
param s15 =1; #Yes
param s16 =0; #No
param s17 =0; #No
param s18 =1; #Yes


#Constraints for non-MEC, non-trivial states
subject to const_s0_a0: s0 >= 0.5 + 0.5*s8;
subject to const_s0_a1: s0 <= 0.5 + 0.5*s8;
subject to const_s2_a0: s2 >= 0.5 + 0;
subject to const_s2_a1: s2 <= 0.5 + 0;
subject to const_s3_a0: s3 >= 0 + 0.5;
subject to const_s3_a1: s3 <= 0 + 0.5;
subject to const_s8_a1: s8 >= s2;
subject to const_s8_a2: s8 >= s3;

#MEC-Constraints

#Objective function
minimize value: (s8 - ( s2)) * (s8 - ( s3)) + 
0;
