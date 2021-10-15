#variable declarations
var s0 >=0, <=1;
param s1 =0; #No
param s2 =0; #No
param s3 =0; #No
param s4 =0; #No
var s5 >=0, <=1;
var s6 >=0, <=1;
var s7 >=0, <=1;
param s8 =0; #No
var s9 >=0, <=1;
var s10 >=0, <=1;
var s11 >=0, <=1;
var s12 >=0, <=1;
param s13 =0; #No
var s14 >=0, <=1;
var s15 >=0, <=1;
var s16 >=0, <=1;
var s17 >=0, <=1;
var s18 >=0, <=1;
var s19 >=0, <=1;
param s20 =0; #No
var s21 >=0, <=1; #MECState
param s22 =0; #No
var s23 >=0, <=1; #MECState
var s24 >=0, <=1; #MECState
param s25 =0; #No
var s26 >=0, <=1; #MECState
var s27 >=0, <=1; #MECState
var s28 >=0, <=1; #MECState
param s29 =0; #No
var s30 >=0, <=1; #MECState
param s31 =0; #No
param s32 =0; #No
var s33 >=0, <=1; #MECState
param s34 =0; #No
param s35 =0; #No
var s36 >=0, <=1; #MECState
var s37 >=0, <=1; #MECState
param s38 =0; #No
param s39 =0; #No
param s40 =0; #No
var s41 >=0, <=1; #MECState
var s42 >=0, <=1; #MECState
param s43 =0; #No
var s44 >=0, <=1; #MECState
param s45 =0; #No
param s46 =0; #No
var s47 >=0, <=1; #MECState
param s48 =0; #No
param s49 =0; #No
var s50 >=0, <=1; #MECState
var s51 >=0, <=1; #MECState
var s52 >=0, <=1; #MECState
var s53 >=0, <=1; #MECState
param s54 =0; #No
param s55 =0; #No
param s56 =0; #No
param s57 =0; #No
param s58 =0; #No
var s59 >=0, <=1; #MECState
var s60 >=0, <=1; #MECState
param s61 =0; #No
param s62 =0; #No
var s63 >=0, <=1; #MECState
var s64 >=0, <=1; #MECState
var s65 >=0, <=1; #MECState
var s66 >=0, <=1; #MECState
var s67 >=0, <=1; #MECState
var s68 >=0, <=1; #MECState
param s69 =0; #No
param s70 =0; #No
param s71 =0; #No
var s72 >=0, <=1; #MECState
var s73 >=0, <=1; #MECState
param s74 =0; #No
var s75 >=0, <=1; #MECState
var s76 >=0, <=1; #MECState
param s77 =0; #No
var s78 >=0, <=1; #MECState
var s79 >=0, <=1; #MECState
param s80 =0; #No
var s81 >=0, <=1; #MECState
var s82 >=0, <=1; #MECState
param s83 =0; #No
param s84 =0; #No
param s85 =0; #No
param s86 =0; #No
var s87 >=0, <=1; #MECState
var s88 >=0, <=1; #MECState
var s89 >=0, <=1; #MECState
var s90 >=0, <=1; #MECState
var s91 >=0, <=1; #MECState
param s92 =0; #No
var s93 >=0, <=1; #MECState
var s94 >=0, <=1; #MECState
var s95 >=0, <=1; #MECState
param s96 =0; #No
var s97 >=0, <=1; #MECState
param s98 =0; #No
var s99 >=0, <=1;
param s100 =0; #No
param s101 =1; #Yes


#Constraints for non-MEC, non-trivial states
subject to const_s0_a1: s0 >= 0.6*s0 + 0;
subject to const_s0_a2: s0 >= 0 + 0;
subject to const_s0_a3: s0 >= 0.4*s0 + 0;
subject to const_s0_a4: s0 >= 0 + 0 + 0;
subject to const_s0_a5: s0 >= 0.3*s0 + 0 + 0.6*s5;
subject to const_s5_a1: s5 <= s26;
subject to const_s5_a2: s5 <= 0 + 0 + 0.6*s27;
subject to const_s5_a3: s5 <= s28;
subject to const_s5_a4: s5 <= 0.2*s21 + 0;
subject to const_s5_a5: s5 <= 0.7*s17 + 0.3*s30;
subject to const_s6_a1: s6 >= 0;
subject to const_s6_a2: s6 >= 0 + 0.1*s18 + 0.8*s6;
subject to const_s6_a3: s6 >= 0.3*s33 + 0.2*s5 + 0.5*s10;
subject to const_s6_a4: s6 >= 0 + 0;
subject to const_s6_a5: s6 >= 0 + 0.1*s19;
subject to const_s7_a1: s7 >= 0.1*s16 + 0.9*s36;
subject to const_s7_a2: s7 >= 0.2*s0 + 0 + 0.4*s37 + 0.2*s7;
subject to const_s7_a3: s7 >= 0 + 0.7*s11;
subject to const_s7_a4: s7 >= 0 + 0;
subject to const_s7_a5: s7 >= 0 + 0.4*s12;
subject to const_s9_a1: s9 >= 0;
subject to const_s9_a2: s9 >= s47;
subject to const_s9_a3: s9 >= 0;
subject to const_s9_a4: s9 >= 0 + 0.3*s11;
subject to const_s9_a5: s9 >= 0.9*s50 + 0.1*s30;
subject to const_s10_a1: s10 >= s51;
subject to const_s10_a2: s10 >= 0.8*s52 + 0.1*s7 + 0.1*s12;
subject to const_s10_a3: s10 >= 0.6*s53 + 0;
subject to const_s10_a4: s10 >= 0 + 0 + 0.2*s23;
subject to const_s10_a5: s10 >= 0 + 0 + 0.1*s27 + 0;
subject to const_s11_a1: s11 >= 0 + 0 + 0;
subject to const_s11_a2: s11 >= 0 + 0;
subject to const_s11_a3: s11 >= 0 + 0;
subject to const_s11_a4: s11 >= 0.4*s33 + 0.3*s10 + 0.3*s59;
subject to const_s11_a5: s11 >= s60;
subject to const_s12_a1: s12 >= 0.3*s36 + 0;
subject to const_s12_a2: s12 >= 0 + 0.1*s5 + 0 + 0.1*s14;
subject to const_s12_a3: s12 >= s63;
subject to const_s12_a4: s12 >= 0.6*s64 + 0 + 0;
subject to const_s12_a5: s12 >= s65;
subject to const_s14_a1: s14 >= 0 + 0 + 0;
subject to const_s14_a2: s14 >= 0.5*s72 + 0.2*s24 + 0;
subject to const_s14_a3: s14 >= 0.1*s19 + 0 + 0.8*s73;
subject to const_s14_a4: s14 >= 0.1*s17 + 0;
subject to const_s14_a5: s14 >= 0.2*s6 + 0.4*s9 + 0.3*s75 + 0;
subject to const_s15_a1: s15 >= 0.2*s42 + 0.8*s76;
subject to const_s15_a2: s15 >= 0.5*s26 + 0;
subject to const_s15_a3: s15 >= 0.8*s18 + 0.1*s37 + 0.1*s78;
subject to const_s15_a4: s15 >= 0.1*s33 + 0.1*s36 + 0 + 0.6*s79;
subject to const_s15_a5: s15 >= 0 + 0.4*s18 + 0.2*s73 + 0.1*s59;
subject to const_s16_a1: s16 >= 0.2*s64 + 0.1*s81 + 0.5*s50 + 0;
subject to const_s16_a2: s16 >= 0.3*s0 + 0.6*s82 + 0;
subject to const_s16_a3: s16 >= 0 + 0 + 0.1*s37 + 0 + 0;
subject to const_s16_a4: s16 >= 0 + 0.1*s42;
subject to const_s16_a5: s16 >= 0 + 0.4*s23;
subject to const_s17_a1: s17 >= 0.1*s66 + 0 + 0.1*s41;
subject to const_s17_a2: s17 >= 0.5*s87 + 0 + 0.2*s72;
subject to const_s17_a3: s17 >= 0.9*s88 + 0;
subject to const_s17_a4: s17 >= 0.9*s89 + 0.1*s63;
subject to const_s17_a5: s17 >= 0 + 0.7*s90 + 0;
subject to const_s18_a1: s18 <= 0 + 0.6*s91 + 0.2*s76;
subject to const_s18_a2: s18 <= 0.9*s75 + 0;
subject to const_s18_a3: s18 <= s93;
subject to const_s18_a4: s18 <= s94;
subject to const_s18_a5: s18 <= s95;
subject to const_s19_a1: s19 <= 0 + 0.8*s30;
subject to const_s19_a2: s19 <= 0.8*s97 + 0;
subject to const_s19_a3: s19 <= 0 + 0 + 0.4*s21 + 0.1*s28 + 0;
subject to const_s19_a4: s19 <= s99;
subject to const_s99_a1: s99 <= 0 + 0.9;
subject to const_s99_a2: s99 <= 0 + 0.01;

#MEC-Constraints
param tmpBestExit0 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s26: s26 = max( tmpBestExit0);
param tmpBestExit1 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s28: s28 = max( tmpBestExit1);
param tmpBestExit2 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s89: s89 = max( tmpBestExit2);
param tmpBestExit3 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s63: s63 = max( tmpBestExit3);
param tmpBestExit4 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s30: s30 = max( tmpBestExit4);
param tmpBestExit5 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s33: s33 = max( tmpBestExit5);
param tmpBestExit6 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s51: s51 = max( tmpBestExit6);
param tmpBestExit7 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s47: s47 = max( tmpBestExit7);
param tmpBestExit8 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s50: s50 = max( tmpBestExit8);
param tmpBestExit9 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s59: s59 = max( tmpBestExit9);
param tmpBestExit10 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s60: s60 = max( tmpBestExit10);
param tmpBestExit11 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s65: s65 = max( tmpBestExit11);
param tmpBestExit12 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s42: s42 = max( tmpBestExit12);
param tmpBestExit13 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s76: s76 = max( tmpBestExit13);
param tmpBestExit14 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s93: s93 = max( tmpBestExit14);
param tmpBestExit15 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s94: s94 = max( tmpBestExit15);
param tmpBestExit16 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s95: s95 = max( tmpBestExit16);
param tmpBestExit17 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s37: s37 = max( tmpBestExit17);
param tmpBestExit18 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s78: s78 = max( tmpBestExit18);
param tmpBestExit19 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s21: s21 = max( tmpBestExit19);
param tmpBestExit20 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s23: s23 = max( tmpBestExit20);
param tmpBestExit21 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s24: s24 = max( tmpBestExit21);
param tmpBestExit22 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s27: s27 = max( tmpBestExit22);
param tmpBestExit23 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s36: s36 = max( tmpBestExit23);
param tmpBestExit24 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s41: s41 = max( tmpBestExit24);
param tmpBestExit25 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s44: s44 = max( tmpBestExit25);
param tmpBestExit26 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s52: s52 = max( tmpBestExit26);
param tmpBestExit27 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s53: s53 = max( tmpBestExit27);
param tmpBestExit28 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s64: s64 = max( tmpBestExit28);
param tmpBestExit29 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s66: s66 = max( tmpBestExit29);
param tmpBestExit30 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s67: s67 = max( tmpBestExit30);
param tmpBestExit31 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s68: s68 = max( tmpBestExit31);
param tmpBestExit32 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s72: s72 = max( tmpBestExit32);
param tmpBestExit33 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s73: s73 = max( tmpBestExit33);
param tmpBestExit34 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s75: s75 = max( tmpBestExit34);
param tmpBestExit35 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s79: s79 = max( tmpBestExit35);
param tmpBestExit36 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s81: s81 = max( tmpBestExit36);
param tmpBestExit37 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s82: s82 = max( tmpBestExit37);
param tmpBestExit38 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s87: s87 = max( tmpBestExit38);
param tmpBestExit39 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s88: s88 = max( tmpBestExit39);
param tmpBestExit40 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s90: s90 = max( tmpBestExit40);
param tmpBestExit41 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s91: s91 = max( tmpBestExit41);
param tmpBestExit42 = 0.99 * s100 +  0.01 * s101;
subject to const_BestExit_s97: s97 = max( tmpBestExit42);

#Objective function
minimize value: (s0 - ( 0.6*s0 +  0)) * (s0 - ( 0.6*s0 +  0)) * (s0 - ( 0 +  0)) * (s0 - ( 0.4*s0 +  0)) * (s0 - ( 0 +  0 +  0)) * (s0 - ( 0.3*s0 +  0 +  0.6*s5)) + 
(s5 - ( s26)) * (s5 - ( s26)) * (s5 - ( 0 +  0 +  0.6*s27)) * (s5 - ( s28)) * (s5 - ( 0.2*s21 +  0)) * (s5 - ( 0.7*s17 +  0.3*s30)) + 
(s6 - ( 0)) * (s6 - ( 0)) * (s6 - ( 0 +  0.1*s18 +  0.8*s6)) * (s6 - ( 0.3*s33 +  0.2*s5 +  0.5*s10)) * (s6 - ( 0 +  0)) * (s6 - ( 0 +  0.1*s19)) + 
(s7 - ( 0.1*s16 +  0.9*s36)) * (s7 - ( 0.1*s16 +  0.9*s36)) * (s7 - ( 0.2*s0 +  0 +  0.4*s37 +  0.2*s7)) * (s7 - ( 0 +  0.7*s11)) * (s7 - ( 0 +  0)) * (s7 - ( 0 +  0.4*s12)) + 
(s9 - ( 0)) * (s9 - ( 0)) * (s9 - ( s47)) * (s9 - ( 0)) * (s9 - ( 0 +  0.3*s11)) * (s9 - ( 0.9*s50 +  0.1*s30)) + 
(s10 - ( s51)) * (s10 - ( s51)) * (s10 - ( 0.8*s52 +  0.1*s7 +  0.1*s12)) * (s10 - ( 0.6*s53 +  0)) * (s10 - ( 0 +  0 +  0.2*s23)) * (s10 - ( 0 +  0 +  0.1*s27 +  0)) + 
(s11 - ( 0 +  0 +  0)) * (s11 - ( 0 +  0 +  0)) * (s11 - ( 0 +  0)) * (s11 - ( 0 +  0)) * (s11 - ( 0.4*s33 +  0.3*s10 +  0.3*s59)) * (s11 - ( s60)) + 
(s12 - ( 0.3*s36 +  0)) * (s12 - ( 0.3*s36 +  0)) * (s12 - ( 0 +  0.1*s5 +  0 +  0.1*s14)) * (s12 - ( s63)) * (s12 - ( 0.6*s64 +  0 +  0)) * (s12 - ( s65)) + 
(s14 - ( 0 +  0 +  0)) * (s14 - ( 0 +  0 +  0)) * (s14 - ( 0.5*s72 +  0.2*s24 +  0)) * (s14 - ( 0.1*s19 +  0 +  0.8*s73)) * (s14 - ( 0.1*s17 +  0)) * (s14 - ( 0.2*s6 +  0.4*s9 +  0.3*s75 +  0)) + 
(s15 - ( 0.2*s42 +  0.8*s76)) * (s15 - ( 0.2*s42 +  0.8*s76)) * (s15 - ( 0.5*s26 +  0)) * (s15 - ( 0.8*s18 +  0.1*s37 +  0.1*s78)) * (s15 - ( 0.1*s33 +  0.1*s36 +  0 +  0.6*s79)) * (s15 - ( 0 +  0.4*s18 +  0.2*s73 +  0.1*s59)) + 
(s16 - ( 0.2*s64 +  0.1*s81 +  0.5*s50 +  0)) * (s16 - ( 0.2*s64 +  0.1*s81 +  0.5*s50 +  0)) * (s16 - ( 0.3*s0 +  0.6*s82 +  0)) * (s16 - ( 0 +  0 +  0.1*s37 +  0 +  0)) * (s16 - ( 0 +  0.1*s42)) * (s16 - ( 0 +  0.4*s23)) + 
(s17 - ( 0.1*s66 +  0 +  0.1*s41)) * (s17 - ( 0.1*s66 +  0 +  0.1*s41)) * (s17 - ( 0.5*s87 +  0 +  0.2*s72)) * (s17 - ( 0.9*s88 +  0)) * (s17 - ( 0.9*s89 +  0.1*s63)) * (s17 - ( 0 +  0.7*s90 +  0)) + 
(s18 - ( 0 +  0.6*s91 +  0.2*s76)) * (s18 - ( 0 +  0.6*s91 +  0.2*s76)) * (s18 - ( 0.9*s75 +  0)) * (s18 - ( s93)) * (s18 - ( s94)) * (s18 - ( s95)) + 
(s19 - ( 0 +  0.8*s30)) * (s19 - ( 0.8*s97 +  0)) * (s19 - ( 0 +  0 +  0.4*s21 +  0.1*s28 +  0)) * (s19 - ( s99)) + 
(s99 - ( 0 +  0.9)) * (s99 - ( 0 +  0.01)) + 
0;
