#variable declarations
var s1 >=0, <=1;
var s4421 >=0, <=1; #MECState


#Constraints for non-MEC, non-trivial states
subject to const_s1_a1: s1 <= 0.001*0.999 + 0.999*s4421;
subject to const_s1_a2: s1 <= 0.99*s4421 + 0.01*0.9999999989989999;

#MEC-Constraints
var tmpBestExit0 >= 0, <= 1;
subject to const_tmpBestExit0_LEQ: tmpBestExit0 <=  0.999 * 1.0 +  0.001 * s1;
subject to const_tmpBestExit0_GEQ: tmpBestExit0 >=  0.999 * 1.0 +  0.001 * s1;
var tmpBestExit1 >= 0, <= 1;
subject to const_tmpBestExit1_LEQ: tmpBestExit1 <=  0.001 * s1 +  0.999 * 1.0;
subject to const_tmpBestExit1_GEQ: tmpBestExit1 >=  0.001 * s1 +  0.999 * 1.0;
var tmpBestExit2 >= 0, <= 1;
subject to const_tmpBestExit2_LEQ: tmpBestExit2 <=  0.001 * s1 +  0.999 * 1.0;
subject to const_tmpBestExit2_GEQ: tmpBestExit2 >=  0.001 * s1 +  0.999 * 1.0;
var tmpBestExit3 >= 0, <= 1;
subject to const_tmpBestExit3_LEQ: tmpBestExit3 <=  0.001 * s1 +  0.999 * 0.9999989990000009;
subject to const_tmpBestExit3_GEQ: tmpBestExit3 >=  0.001 * s1 +  0.999 * 0.9999989990000009;
var tmpBestExit4 >= 0, <= 1;
subject to const_tmpBestExit4_LEQ: tmpBestExit4 <=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_tmpBestExit4_GEQ: tmpBestExit4 >=  0.01 * s1 +  0.99 * 0.9999989990000009;
var tmpBestExit5 >= 0, <= 1;
subject to const_tmpBestExit5_LEQ: tmpBestExit5 <=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_tmpBestExit5_GEQ: tmpBestExit5 >=  0.01 * s1 +  0.99 * 0.9999989990000009;
var tmpBestExit6 >= 0, <= 1;
subject to const_tmpBestExit6_LEQ: tmpBestExit6 <=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_tmpBestExit6_GEQ: tmpBestExit6 >=  0.01 * s1 +  0.99 * 0.9999989990000009;
var tmpBestExit7 >= 0, <= 1;
subject to const_tmpBestExit7_LEQ: tmpBestExit7 <=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_tmpBestExit7_GEQ: tmpBestExit7 >=  0.01 * s1 +  0.99 * 0.9999989990000009;
var tmpBestExit8 >= 0, <= 1;
subject to const_tmpBestExit8_LEQ: tmpBestExit8 <=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_tmpBestExit8_GEQ: tmpBestExit8 >=  0.01 * s1 +  0.99 * 0.9999989990000009;
var tmpBestExit9 >= 0, <= 1;
subject to const_tmpBestExit9_LEQ: tmpBestExit9 <=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_tmpBestExit9_GEQ: tmpBestExit9 >=  0.01 * s1 +  0.99 * 0.9999989990000009;
subject to const_BestExit_s4421_Min: s4421 >= max( tmpBestExit0, tmpBestExit1, tmpBestExit2, tmpBestExit3, tmpBestExit4, tmpBestExit5, tmpBestExit6, tmpBestExit7, tmpBestExit8, tmpBestExit9);
subject to const_BestExit_s4421_Max: s4421 <= max( tmpBestExit0, tmpBestExit1, tmpBestExit2, tmpBestExit3, tmpBestExit4, tmpBestExit5, tmpBestExit6, tmpBestExit7, tmpBestExit8, tmpBestExit9);

#Objective function
minimize value: (s1 - ( 0.001*0.999 +  0.999*s4421)) * (s1 - ( 0.99*s4421 +  0.01*0.9999999989989999)) + 
0;
