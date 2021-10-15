#variable declarations
var s1 >=0, <=1;
var s17477 >=0, <=1; #MECState


#Constraints for non-MEC, non-trivial states
subject to const_s1_a1: s1 <= 0.001*0.99999899000101 + 0.999*s17477;
subject to const_s1_a2: s1 <= 0.99*s17477 + 0.01*0.99;

#MEC-Constraints
var tmpBestExit0 >= 0, <= 1;
subject to const_tmpBestExit0_LEQ: tmpBestExit0 <=  0.001 * s1 +  0.999 * 0.9998990001009999;
subject to const_tmpBestExit0_GEQ: tmpBestExit0 >=  0.001 * s1 +  0.999 * 0.9998990001009999;
var tmpBestExit1 >= 0, <= 1;
subject to const_tmpBestExit1_LEQ: tmpBestExit1 <=  0.001 * s1 +  0.999 * 0.9998990001009999;
subject to const_tmpBestExit1_GEQ: tmpBestExit1 >=  0.001 * s1 +  0.999 * 0.9998990001009999;
var tmpBestExit2 >= 0, <= 1;
subject to const_tmpBestExit2_LEQ: tmpBestExit2 <=  0.001 * s1 +  0.999 * 0.9998990001009999;
subject to const_tmpBestExit2_GEQ: tmpBestExit2 >=  0.001 * s1 +  0.999 * 0.9998990001009999;
var tmpBestExit3 >= 0, <= 1;
subject to const_tmpBestExit3_LEQ: tmpBestExit3 <=  0.001 * s1 +  0.999 * 0.9998990001009999;
subject to const_tmpBestExit3_GEQ: tmpBestExit3 >=  0.001 * s1 +  0.999 * 0.9998990001009999;
var tmpBestExit4 >= 0, <= 1;
subject to const_tmpBestExit4_LEQ: tmpBestExit4 <=  0.01 * s1 +  0.99 * 0.9998990001009999;
subject to const_tmpBestExit4_GEQ: tmpBestExit4 >=  0.01 * s1 +  0.99 * 0.9998990001009999;
var tmpBestExit5 >= 0, <= 1;
subject to const_tmpBestExit5_LEQ: tmpBestExit5 <=  0.01 * s1 +  0.99 * 0.9998990001009999;
subject to const_tmpBestExit5_GEQ: tmpBestExit5 >=  0.01 * s1 +  0.99 * 0.9998990001009999;
var tmpBestExit6 >= 0, <= 1;
subject to const_tmpBestExit6_LEQ: tmpBestExit6 <=  0.01 * s1 +  0.99 * 0.9998990001009999;
subject to const_tmpBestExit6_GEQ: tmpBestExit6 >=  0.01 * s1 +  0.99 * 0.9998990001009999;
var tmpBestExit7 >= 0, <= 1;
subject to const_tmpBestExit7_LEQ: tmpBestExit7 <=  0.01 * s1 +  0.99 * 0.9998990001009999;
subject to const_tmpBestExit7_GEQ: tmpBestExit7 >=  0.01 * s1 +  0.99 * 0.9998990001009999;
var tmpBestExit8 >= 0, <= 1;
subject to const_tmpBestExit8_LEQ: tmpBestExit8 <=  0.01 * s1 +  0.99 * 1.0;
subject to const_tmpBestExit8_GEQ: tmpBestExit8 >=  0.01 * s1 +  0.99 * 1.0;
var tmpBestExit9 >= 0, <= 1;
subject to const_tmpBestExit9_LEQ: tmpBestExit9 <=  0.01 * s1 +  0.99 * 1.0;
subject to const_tmpBestExit9_GEQ: tmpBestExit9 >=  0.01 * s1 +  0.99 * 1.0;
var tmpBestExit10 >= 0, <= 1;
subject to const_tmpBestExit10_LEQ: tmpBestExit10 <=  0.01 * s1 +  0.99 * 1.0;
subject to const_tmpBestExit10_GEQ: tmpBestExit10 >=  0.01 * s1 +  0.99 * 1.0;
var tmpBestExit11 >= 0, <= 1;
subject to const_tmpBestExit11_LEQ: tmpBestExit11 <=  0.01 * s1 +  0.99 * 0.9998990001009999;
subject to const_tmpBestExit11_GEQ: tmpBestExit11 >=  0.01 * s1 +  0.99 * 0.9998990001009999;
subject to const_BestExit_s17477_Min: s17477 >= max( tmpBestExit0, tmpBestExit1, tmpBestExit2, tmpBestExit3, tmpBestExit4, tmpBestExit5, tmpBestExit6, tmpBestExit7, tmpBestExit8, tmpBestExit9, tmpBestExit10, tmpBestExit11);
subject to const_BestExit_s17477_Max: s17477 <= max( tmpBestExit0, tmpBestExit1, tmpBestExit2, tmpBestExit3, tmpBestExit4, tmpBestExit5, tmpBestExit6, tmpBestExit7, tmpBestExit8, tmpBestExit9, tmpBestExit10, tmpBestExit11);

#Objective function
minimize value: (s1 - ( 0.001*0.99999899000101 +  0.999*s17477)) * (s1 - ( 0.99*s17477 +  0.01*0.99)) + 
0;
