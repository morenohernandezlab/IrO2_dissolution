 function products = water_radiolysis(t,reactants,generationDueToGValues,k)

% Water-Cl-Ce Radiolysis Model
% Modified by Ivan A. Moreno-Hernandez to include Cl and Ce reactions, 2020
% Modified by Ivan A. Moreno-Hernandez and S. Avery Vigil to include Fe reactions, 2023
% Modified by Ivan A. Moreno-Hernandez and Matteo Fratarcangeli to expand the Fe-Cl(ox) interactions, 2024

%{
    The MIT License (MIT)
    
    Copyright (c) 2014 Brian J. Mendel, Nicholas M. Schneider
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
        
        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
        
%}



% This function, when used in conjunction with an ODE solver (ode23s and
% ode15s work best).
%
% Concentrations are in MICROMOLAR, and time is in SECONDS.
% 
% 
%"t" is the time span [tinitial tfinal] in SECONDS, "reactants" is a row  
% vector of 16 elements specifying the initial concentrations of radiolysis
% products in MICROMOLAR.



%Molar to micromolar conversion factor%
molarToMicromolar = 10^6;



%%convert from molar to micromolar
k = k./(molarToMicromolar);






%Reactant concentrations (micromolar)%
electrons = reactants(1); 
H_ions    = reactants(2);
OH_ions   = reactants(3);
H2O2      = reactants(4);
HO2_ions  = reactants(5);
H_mono    = reactants(6);
OH        = reactants(7);
O_ions    = reactants(8);
HO2       = reactants(9);
O2_ions   = reactants(10);
O2        = reactants(11);
H2        = reactants(12);
O3_ions   = reactants(13);
O3        = reactants(14);
HO3       = reactants(15);
H2O       = reactants(16);
Cl_ions   = reactants(17);
ClOH_ions = reactants(18);
Cl_mono   = reactants(19);
Cl2_ions  = reactants(20);
Cl3_ions  = reactants(21);
Cl2       = reactants(22);
Fe2_ions  = reactants(23);
Fe3_ions  = reactants(24);
Fe4_ions  = reactants(25);
Fe5       = reactants(26);
Fe6_ions  = reactants(27);
Fe3_OH_ions =  reactants(28);
Fe3_OH2_ions =  reactants(29);
Fe3_dimer_ions = reactants(30);
Fe2_OH_ions    = reactants(31);
Fe3_HO2_ions = reactants(32);
Fe3_H2O3_ions = reactants(33);
O2_2m_ions = reactants(34);
Fe2_Cl_ions = reactants(35);
Fe3_Cl_ions = reactants(36);
Fe3_Cl2_ions = reactants(37);

%Reaction Set%
% Note: H2O is divided out as a reactant
r = zeros(73,1);
r(1)  = k(1)*H_ions*OH_ions;
r(2)  = k(2);
r(3)  = k(3)*H2O2;
r(4)  = k(4)*H_ions*HO2_ions;
r(5)  = k(5)*H2O2*OH_ions;
r(6)  = k(6)*HO2_ions;
r(7)  = k(7)*electrons*H2O;
r(8)  = k(8)*H_mono*OH_ions;
r(9)  = k(9)*H_mono;
r(10) = k(10)*electrons*H_ions;
r(11) = k(11)*OH*OH_ions;
r(12) = k(12)*O_ions;
r(13) = k(13)*OH;
r(14) = k(14)*O_ions*H_ions;
r(15) = k(15)*HO2;
r(16) = k(16)*O2_ions*H_ions;
r(17) = k(17)*HO2*OH_ions;
r(18) = k(18)*O2_ions;
r(19) = k(19)*electrons*OH;
r(20) = k(20)*electrons*H2O2;
r(21) = k(21)*electrons*O2_ions;
r(22) = k(22)*electrons*HO2;
r(23) = k(23)*electrons*O2;
r(24) = k(24)*electrons^2;
r(25) = k(25)*electrons*H_mono;
r(26) = k(26)*electrons*HO2_ions;
r(27) = k(27)*electrons*O_ions;
r(28) = k(28)*electrons*O3_ions;
r(29) = k(29)*electrons*O3;
r(30) = k(30)*H_mono*H2O; 
r(31) = k(31)*H_mono*O_ions;
r(32) = k(32)*H_mono*HO2_ions;
r(33) = k(33)*H_mono*O3_ions;
r(34) = k(34)*H_mono^2;
r(35) = k(35)*H_mono*OH;
r(36) = k(36)*H_mono*H2O2;
r(37) = k(37)*H_mono*O2;
r(38) = k(38)*H_mono*HO2;
r(39) = k(39)*H_mono*O2_ions;
r(40) = k(40)*H_mono*O3;
r(41) = k(41)*OH^2;
r(42) = k(42)*OH*HO2; 
r(43) = k(43)*OH*O2_ions;
r(44) = k(44)*OH*H2;
r(45) = k(45)*OH*H2O2;
r(46) = k(46)*OH*O_ions;
r(47) = k(47)*OH*HO2_ions;
r(48) = k(48)*OH*O3_ions;
r(49) = k(49)*OH*O3_ions;
r(50) = k(50)*OH*O3;
r(51) = k(51)*HO2*O2_ions;
r(52) = k(52)*HO2^2;
r(53) = k(53)*HO2*O_ions;
r(54) = k(54)*HO2*H2O2;
r(55) = k(55)*HO2*HO2_ions;
r(56) = k(56)*HO2*O3_ions;
r(57) = k(57)*HO2*O3;
r(58) = k(58)*O2_ions^2;
r(59) = k(59)*O2_ions*O_ions;
r(60) = k(60)*O2_ions*H2O2;
r(61) = k(61)*O2_ions*HO2_ions;
r(62) = k(62)*O2_ions*O3_ions;
r(63) = k(63)*O2_ions*O3;
r(64) = k(64)*O_ions^2;
r(65) = k(65)*O_ions*O2;
r(66) = k(66)*O_ions*H2;
r(67) = k(67)*O_ions*H2O2;
r(68) = k(68)*O_ions*HO2_ions;
r(69) = k(69)*O_ions*O3_ions;
r(70) = k(70)*O_ions*O3;
r(71) = k(71)*O3_ions;
r(72) = k(72)*O3_ions*H_ions;
r(73) = k(73)*HO3;
r(74) = k(74)*OH*Cl_ions;
r(75) = k(75)*ClOH_ions;
r(76) = k(76)*Cl_mono*Cl_ions;
r(77) = k(77)*H_ions*ClOH_ions;
r(78) = k(78)*ClOH_ions;
r(79) = k(79)*Cl_mono*OH;
r(80) = k(80)*Cl2_ions;
r(81) = k(81)*Cl2_ions*Cl2_ions;
r(82) = k(82)*Cl_mono*Cl2_ions;
r(83) = k(83)*Cl_ions*Cl2;
r(84) = k(84)*Cl3_ions;
r(85) = k(85)*Cl_mono*Cl_mono;
r(86) = k(86)*electrons*Cl_mono;
r(87) = k(87)*electrons*Cl2_ions;
r(88) = k(88)*electrons*Cl3_ions;
r(89) = k(89)*H_mono*Cl_mono;
r(90) = k(90)*H_mono*Cl2_ions;
r(91) = k(91)*H_mono*Cl3_ions;
r(92) = k(92)*HO2*Cl2_ions;
r(93) = k(93)*HO2*Cl2;
r(94) = k(94)*Fe2_ions*OH;
r(95) = k(95)*Fe3_ions*OH;
r(96) = k(96)*Fe4_ions*OH;
r(97) = k(97)*Fe5*OH;
r(98) = k(98)*Fe3_ions*H_mono;
r(99) = k(99)*Fe4_ions*H_mono;
r(100) = k(100)*Fe5*H_mono;
r(101) = k(101)*Fe6_ions*H_mono;
r(102) = k(102)*Fe3_ions*electrons;
r(103) = k(103)*Fe4_ions*electrons;
r(104) = k(104)*Fe5*electrons;
r(105) = k(105)*Fe6_ions*electrons;
r(106) = k(106)*Fe6_ions*Fe6_ions; %Self decay of Fe6+
r(107) = k(107)*Fe5*Fe5; %Self-decay of Fe5+
r(108) = k(108)*Fe5*H2O2;
r(109) = k(109)*Fe4_ions*Fe4_ions;
r(110) = k(110)*Fe4_ions*H2O2;
r(111) = k(111)*Fe3_ions;
r(112) = k(112)*H_ions*Fe3_OH_ions;
r(113) = k(113)*Fe3_ions;
r(114) = k(114)*Fe3_OH2_ions*H_ions*H_ions;
r(115) = k(115)*Fe3_ions*Fe3_ions;
r(116) = k(116)*Fe3_dimer_ions*H_ions*H_ions;
r(117) = k(117)*Fe2_ions;
r(118) = k(118)*Fe2_OH_ions*H_ions;
r(119) = k(119)*Fe2_ions*H2O2;
r(120) = k(120)*Fe2_OH_ions*H2O2;
r(121) = k(121)*Fe3_ions*H2O2;
r(122) = k(122)*Fe3_HO2_ions*H_ions;
r(123) = k(123)*Fe3_OH_ions*H2O2;
r(124) = k(124)*Fe3_H2O3_ions*H_ions;
r(125) = k(125)*Fe3_HO2_ions;
r(126) = k(126)*Fe3_H2O3_ions;
r(127) = k(127)*Fe2_OH_ions*OH;
r(128) = k(128)*Fe2_ions*HO2;
r(129) = k(129)*Fe2_ions*O2_ions;
r(130) = k(130)*O2_2m_ions*H_ions;
r(131) = k(131)*Fe3_ions*HO2;
r(132) = k(132)*Fe3_ions*O2_ions;
r(133) = k(133)*Fe2_ions*Cl_ions;
r(134) = k(134)*Fe2_Cl_ions;
r(135) = k(135)*Fe3_ions*Cl_ions;
r(136) = k(136)*Fe3_Cl_ions;
r(137) = k(137)*Fe3_ions*Cl_ions*Cl_ions;
r(138) = k(138)*Fe3_Cl2_ions;
r(139) = k(139)*Fe2_Cl_ions*H2O2;
r(140) = k(140)*Fe2_ions*Cl_mono;
r(141) = k(141)*Fe2_ions*Cl2_ions;
r(142) = k(142)*Fe2_Cl_ions*HO2;
r(143) = k(143)*Fe2_Cl_ions*O2_ions;
r(144) = k(144)*Fe3_Cl_ions*HO2;
r(145) = k(145)*Fe3_Cl2_ions*HO2;
r(146) = k(146)*Fe3_Cl_ions*O2_ions;
r(147) = k(147)*Fe3_Cl2_ions*O2_ions;
r(148) = k(148)*Fe2_Cl_ions*OH;
r(149) = k(149)*Fe3_Cl_ions*OH;
r(150) = k(150)*Fe3_Cl2_ions*OH;
r(151) = k(151)*Fe3_Cl_ions*H_mono;
r(152) = k(152)*Fe3_Cl2_ions*H_mono;
r(153) = k(153)*Fe3_Cl_ions*electrons;
r(154) = k(154)*Fe3_Cl2_ions*electrons;
r(155) = k(155)*Fe2_Cl_ions*Cl_mono;
r(156) = k(156)*Fe2_Cl_ions*Cl2_ions;
r(157) = k(157)*Fe2_ions*Cl3_ions;
r(158) = k(158)*Fe2_ions*Cl2;
r(159) = k(159)*Fe2_Cl_ions*Cl3_ions;
r(160) = k(160)*Fe2_Cl_ions*Cl2;

products = zeros(length(reactants),1);


%%electrons%%
products(1) = -r(7) + r(8) + r(9) - r(10) - r(19) - r(20) - r(21) - r(22)...
    -r(23) - 2*r(24) - r(25) - r(26) - r(27) - r(28) - r(29)-r(86)-r(87)-r(88)...
    -r(102)-r(103)-r(104)-r(105)-r(153)-r(154);

%%H+%%
products(2) = -r(1) + r(2) + r(3) - r(4) + r(9) - r(10) + r(13) - r(14)...
    + r(15) - r(16) + r(49) - r(72)-r(77)+r(89)+r(90)+r(91)+r(92)+r(93)+5*r(95)...
    +r(98)-4*r(99)+r(100)+r(101)-5*r(103)-2*r(106)-3*r(107)-3*r(108)-8*r(109)-3*r(110)+r(111)-r(112)+2*r(113)-2*r(114)+2*r(115)-2*r(116)+r(117)-r(118)+r(121)-r(122)+r(123)-r(124)-r(130)+r(131)+r(144)+r(145)+5*r(149)+5*r(150)+r(151)+r(152);

%%OH-%%
products(3) = -r(1) + r(2) - r(5) + r(6) + r(7) - r(8) - r(11) + r(12) -...
    r(17)+ r(18) + r(19) + r(20) + r(21) + 2*r(24) + r(25) + r(26) + ...
    2*r(27) + 2*r(28) + r(31) + r(32) + r(33) + r(43) + r(47) + r(48) + ...
    r(53) + r(55) + r(56) + 2*r(58) + 2*r(59) + r(60) + r(61) + 2*r(62)...
    + r(64) + r(66) + r(68)+r(78)-r(79)+r(94)+r(95)+r(96)+r(97)+r(119)+2*r(120)+r(126)+2*r(127)+r(139)+r(148)+r(149)+r(150);

%%H2O2%%
products(4) = - r(3) + r(4) -r(5) + r(6) - r(20) - r(36) + r(38) + r(41)...
    - r(45) + r(52) - r(54) + r(58) - r(60) - r(67)+0.3*r(106)+r(107)-r(108)+r(109)-r(110)-r(119)-r(120)-r(121)+r(122)-r(123)+r(124)-r(139);

%%HO2-%%
products(5) = r(3) - r(4) + r(5) - r(6) + r(21) + r(22) - r(26) - r(32)...
    + r(39) + r(46) - r(47) + r(51) - r(55) - r(61) + r(64) - r(68)+r(128)+r(130)+r(142);

%%H%%
products(6) = r(7) - r(8) - r(9) + r(10) - r(25) - r(30) - r(31) - r(32)...
    - r(33) - 2*r(34) - r(35) - r(36) - r(37) - r(38) - r(39) - r(40) + r(44)...
    + r(66)-r(89)-r(90)-r(91)-r(98)-r(99)-r(100)-r(101)-r(151)-r(152);

%%OH%%
products(7) = -r(11) + r(12) - r(13) + r(14) - r(19) + r(20) + r(30) ...
    + r(32) - r(35) + r(36) - 2*r(41) - r(42) - r(43) - r(44) - r(45) -...
    r(46) - r(47) - r(48) - r(49) - r(50) + r(54) + r(55) + r(60) + r(72)...
    + r(73)-r(74)+r(75)-r(94)-r(95)-r(96)-r(97)+r(119)+r(120)-r(127)+r(139)-r(148)-r(149)-r(150);

%%O-%%
products(8) = r(11) - r(12) + r(13) - r(14) + r(26) - r(27) - r(31) - r(46)...
    - r(53) - r(59) + r(61) - 2*r(64) - r(65) - r(66) - r(67) - r(68) -...
    r(69) - r(70) + r(71) ;

%%HO2%%
products(9) = -r(15) + r(16) - r(17) + r(18) - r(22) + r(37) - r(38) - ...
    r(42) + r(45) + r(47) + r(50) - r(51) - 2*r(52) - r(53) - r(54) - r(55)...
    - r(56) - r(57)-r(92)-r(93)+r(125)+r(126)-r(128)-r(131)-r(142)-r(144)-r(145);

%%O2-%%
products(10) = r(15) - r(16) + r(17) - r(18) - r(21) + r(23) - r(39) - ...
    r(43) + 2*r(49) - r(51) - 2*r(58) - r(59) - r(60) - r(61) - r(62) - r(63)...
    + r(67) + r(68) + 2*r(69) + r(70)-r(129)-r(132)-r(143)-r(146)-r(147);

%%O2%%
products(11) = -r(23) + r(28) + r(33) - r(37) + r(42) + r(43) + r(50) ...
    + r(51) + r(52) + r(53) + r(54) + r(55) + 2*r(56) + r(57) + r(58) + r(59)...
    + r(60) + r(61) + 2*r(62) + r(63) - r(65) + r(70) + r(71) + r(72) + r(73)+r(92)+r(93)+0.6*r(106)+r(108)+r(110)+r(131)+r(132)+r(144)+r(145)+r(146)+r(147);

%%H2%%
products(12) = r(24) + r(25) + r(30) + r(34) - r(44) - r(66);

%%O3-%%
products(13) = -r(28) + r(29) - r(33) - r(48) - r(49) - r(56) - r(62) ...
    + r(63) + r(65) - r(69) - r(71) - r(72);

%%O3%%
products(14) = -r(29) - r(40) + r(48) - r(50) - r(57) - r(63) - r(70);

%%HO3%%
products(15) = r(40) + r(57) - r(73);

%%H2O%%
products(16) = r(1) - r(2) + r(5) - r(6) - r(7) + r(8) + r(11) - r(12) + ...
     r(17) - r(18) - r(21) - 2*r(24) - r(25) - r(27) - r(28) - r(30) + r(35)...
     + r(36) + r(42) + r(44) + r(45) + r(54) - 2*r(58) - r(59) - r(62) - r(64) ...
     + r(67)+r(77)-4*r(95)+4*r(99)+4*r(103)+2.2*r(106)+2*r(107)+4*r(108)+6*r(109)+4*r(110)-r(111)+r(112)-2*r(113)+2*r(114)-2*r(115)+2*r(116)-r(117)+r(118)-4*r(149)-4*r(150);
%%Cl-%%
products(17) =-r(74)+r(75)-r(76)+r(80)+r(81)-r(83)+r(84)+r(86)+2*r(87)+r(88)+r(89)+2*r(90)+r(91)+2*r(92)-r(133)+r(134)-r(135)+r(136)-2*r(137)+2*r(138)+r(139)+r(140)+2*r(141)+r(142)+r(143)+r(144)+2*r(145)+r(146)+2*r(147)+r(148)+r(149)+2*r(150)+r(151)+2*r(152)+r(153)+2*r(154)+r(155)+2*r(156)+3*r(157)+2*r(158)+3*r(159)+2*r(160);

%%ClOH-%%
products(18) =r(74)-r(75)-r(77)-r(78)+r(79);

%%Cl%%
products(19) =-r(76)+r(77)+r(78)-r(79)+r(80)-r(82)-2*r(85)-r(86)-r(89)-r(140)-r(155);
 
%%Cl2-%%
products(20) =r(76)-r(80)-2*r(81)-r(82)-r(87)+r(88)-r(90)+r(91)-r(92)+r(93)-r(141)-r(156);

%%Cl3-%%
products(21) =r(81)+r(82)+r(83)-r(84)-r(88)-r(91)-r(157)-r(159);

%%Cl2%%
products(22) =-r(83)+r(84)+r(85)-r(93)-r(158)-r(160);

%%Fe2_ions%%
products(23) =-r(94)+r(98)+r(102)+r(110)-r(117)+r(118)-r(119)+r(125)+r(126)-r(128)-r(129)+r(131)+r(132)-r(133)+r(134)-r(140)-r(141)+r(144)+r(145)+r(146)+r(147)+r(151)+r(152)+r(153)+r(154)-2*r(157)-2*r(158);

%%Fe3_ions%%
products(24) =r(94)-r(95)-r(98)+r(99)-r(102)+r(103)+r(106)+r(107)+r(108)+2*r(109)-r(111)+r(112)-r(113)+r(114)-2*r(115)+2*r(116)+r(119)+r(120)-r(121)+r(122)+r(127)+r(128)+r(129)-r(131)-r(132)-r(135)+r(136)-r(137)+r(138)+r(139)+r(140)+r(141)+r(142)+r(143)+r(148)+2*r(157)+2*r(158);

%%Fe4_ions%%
products(25) =r(95)-r(96)-r(99)+r(100)-r(103)+r(104)-2*r(109)-r(110)+r(149)+r(150);

%%Fe5%%
products(26) =r(96)-r(97)-r(100)+r(101)-r(104)+r(105)-r(107)-r(108);

%%Fe6_ions%%
products(27) =r(97)-r(101)-r(105)-r(106);

%%Fe3_OH_ions%%
products(28) =r(111)-r(112)-r(123)+r(124);


%%Fe3_OH2_ions%%
products(29) =r(113)-r(114);

%%Fe3_dimer_ions%%
products(30) =r(115)-r(116);

%%Fe2_OH_ions%%
products(31) =r(117)-r(118)-r(120)-r(127);

%%Fe3_HO2_ions%%
products(32) =r(121)-r(122)-r(125);

%%Fe3_H2O3_ions%%
products(33) =r(123)-r(124)-r(126);


%%O2_2m_ions%%
products(34) =r(129)-r(130)+r(143);


%%Fe2_Cl_ions%%
products(35) =r(133)-r(134)-r(139)-r(142)-r(143)-r(148)-r(155)-r(156)-2*r(159)-2*r(160);

%%Fe3_Cl_ions%%
products(36) =r(135)-r(136)-r(144)-r(146)-r(149)-r(151)-r(153)+r(155)+r(156)+2*r(159)+2*r(160);

%%Fe3_Cl2_ions%%
products(37) =r(137)-r(138)-r(145)-r(147)-r(150)-r(152)-r(154);

products = products + generationDueToGValues;

