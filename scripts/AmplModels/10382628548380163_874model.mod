#variable declarations
var s0 >=0, <=1;
param s1 =1; #Yes
var s2 >=0, <=1;
var s3 >=0, <=1;
param s4 =1; #Yes
var s5 >=0, <=1;
param s6 =0; #No
var s7 >=0, <=1;
param s8 =0; #No
var s9 >=0, <=1;
param s10 =1; #Yes
var s11 >=0, <=1;
var s12 >=0, <=1;
param s13 =1; #Yes
var s14 >=0, <=1;
param s15 =0; #No
var s16 >=0, <=1;
param s17 =0; #No
param s18 =0; #No
var s19 >=0, <=1;
param s20 =0; #No
param s21 =0; #No
var s22 >=0, <=1;
param s23 =0; #No
param s24 =0; #No
param s25 =0; #No
param s26 =0; #No
param s27 =1; #Yes
param s28 =1; #Yes
param s29 =0; #No
param s30 =0; #No
param s31 =0; #No
param s32 =0; #No
param s33 =0; #No
param s34 =0; #No
param s35 =1; #Yes
param s36 =0; #No
param s37 =1; #Yes
param s38 =0; #No
param s39 =0; #No
param s40 =0; #No
param s41 =0; #No
param s42 =0; #No
param s43 =1; #Yes
param s44 =0; #No
param s45 =0; #No
param s46 =0; #No
param s47 =1; #Yes
param s48 =0; #No
param s49 =0; #No
param s50 =0; #No
var s51 >=0, <=1;
param s52 =1; #Yes
var s53 >=0, <=1;
var s54 >=0, <=1;
param s55 =1; #Yes
var s56 >=0, <=1;
param s57 =0; #No
var s58 >=0, <=1;
param s59 =0; #No
var s60 >=0, <=1;
param s61 =1; #Yes
var s62 >=0, <=1;
var s63 >=0, <=1;
param s64 =1; #Yes
var s65 >=0, <=1;
param s66 =0; #No
var s67 >=0, <=1;
param s68 =0; #No
param s69 =1; #Yes
param s70 =1; #Yes
param s71 =0; #No
param s72 =0; #No
param s73 =0; #No
param s74 =0; #No
param s75 =0; #No
param s76 =0; #No
param s77 =1; #Yes
param s78 =0; #No
param s79 =1; #Yes
param s80 =0; #No
param s81 =0; #No
param s82 =0; #No
param s83 =0; #No
param s84 =0; #No
param s85 =1; #Yes
param s86 =0; #No
param s87 =0; #No
param s88 =0; #No
param s89 =1; #Yes
param s90 =0; #No
param s91 =0; #No
param s92 =0; #No
var s93 >=0, <=1;
var s94 >=0, <=1;
param s95 =0; #No
var s96 >=0, <=1;
var s97 >=0, <=1;
param s98 =0; #No
param s99 =0; #No
param s100 =0; #No
param s101 =0; #No


#Constraints for non-MEC, non-trivial states
subject to const_s0_a1: s0 >= s51;
subject to const_s0_a2: s0 >= s60;
subject to const_s0_a3: s0 >= s93;
subject to const_s2_a1: s2 >= s53;
subject to const_s2_a2: s2 >= s62;
subject to const_s3_a1: s3 >= s54;
subject to const_s3_a2: s3 >= s94;
subject to const_s5_a0: s5 >= s56;
subject to const_s5_a1: s5 <= s56;
subject to const_s7_a0: s7 >= s58;
subject to const_s7_a1: s7 <= s58;
subject to const_s9_a1: s9 >= s63;
subject to const_s9_a2: s9 >= s96;
subject to const_s11_a0: s11 >= s65;
subject to const_s11_a1: s11 <= s65;
subject to const_s12_a0: s12 >= s97;
subject to const_s12_a1: s12 <= s97;
subject to const_s14_a0: s14 >= 0 + 0.3333333333333333 + 0;
subject to const_s14_a1: s14 <= 0 + 0.3333333333333333 + 0;
subject to const_s16_a0: s16 >= 0.3333333333333333 + 0 + 0;
subject to const_s16_a1: s16 <= 0.3333333333333333 + 0 + 0;
subject to const_s19_a0: s19 >= s67;
subject to const_s19_a1: s19 <= s67;
subject to const_s22_a0: s22 >= 0 + 0 + 0.3333333333333333;
subject to const_s22_a1: s22 <= 0 + 0 + 0.3333333333333333;
subject to const_s51_a1: s51 >= s9;
subject to const_s51_a2: s51 >= 0;
subject to const_s53_a1: s53 >= s11;
subject to const_s53_a2: s53 >= 0;
subject to const_s54_a1: s54 >= s12;
subject to const_s54_a2: s54 >= 0;
subject to const_s56_a1: s56 >= s14;
subject to const_s56_a2: s56 >= 0;
subject to const_s58_a1: s58 >= s16;
subject to const_s58_a2: s58 >= 0;
subject to const_s60_a1: s60 >= s3;
subject to const_s60_a2: s60 >= 0;
subject to const_s62_a1: s62 >= s5;
subject to const_s62_a2: s62 >= 0;
subject to const_s63_a1: s63 >= s12;
subject to const_s63_a2: s63 >= 0;
subject to const_s65_a1: s65 >= s14;
subject to const_s65_a2: s65 >= 0;
subject to const_s67_a1: s67 >= s22;
subject to const_s67_a2: s67 >= 0;
subject to const_s93_a1: s93 <= 1.0;
subject to const_s93_a2: s93 <= s2;
subject to const_s94_a1: s94 <= 1.0;
subject to const_s94_a2: s94 <= s5;
subject to const_s96_a1: s96 <= 1.0;
subject to const_s96_a2: s96 <= s11;
subject to const_s97_a1: s97 <= 1.0;
subject to const_s97_a2: s97 <= s14;

#MEC-Constraints

#Objective function
minimize value: (s0 - ( s51)) * (s0 - ( s51)) * (s0 - ( s60)) * (s0 - ( s93)) + 
(s2 - ( s53)) * (s2 - ( s62)) + 
(s3 - ( s54)) * (s3 - ( s94)) + 
(s9 - ( s63)) * (s9 - ( s96)) + 
(s51 - ( s9)) * (s51 - ( 0)) + 
(s53 - ( s11)) * (s53 - ( 0)) + 
(s54 - ( s12)) * (s54 - ( 0)) + 
(s56 - ( s14)) * (s56 - ( 0)) + 
(s58 - ( s16)) * (s58 - ( 0)) + 
(s60 - ( s3)) * (s60 - ( 0)) + 
(s62 - ( s5)) * (s62 - ( 0)) + 
(s63 - ( s12)) * (s63 - ( 0)) + 
(s65 - ( s14)) * (s65 - ( 0)) + 
(s67 - ( s22)) * (s67 - ( 0)) + 
(s93 - ( 1.0)) * (s93 - ( s2)) + 
(s94 - ( 1.0)) * (s94 - ( s5)) + 
(s96 - ( 1.0)) * (s96 - ( s11)) + 
(s97 - ( 1.0)) * (s97 - ( s14)) + 
0;
