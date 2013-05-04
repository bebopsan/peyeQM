Point(1) = {0.000000, 0.000000, 0, 1.0};
Point(2) = {1.000000, 0.000000, 0, 1.0};
Point(3) = {1.000000, 1.000000, 0, 1.0};
Point(4) = {0.000000, 1.000000, 0, 1.0};
Point(6) = {2.000000, 0.000000, 0, 1.0};
Point(7) = {2.000000, 1.000000, 0, 1.0};
Point(10) = {3.000000, 0.000000, 0, 1.0};
Point(11) = {3.000000, 1.000000, 0, 1.0};
Point(15) = {1.000000, 2.000000, 0, 1.0};
Point(16) = {0.000000, 2.000000, 0, 1.0};
Point(19) = {2.000000, 2.000000, 0, 1.0};
Point(23) = {3.000000, 2.000000, 0, 1.0};
Point(27) = {1.000000, 3.000000, 0, 1.0};
Point(28) = {0.000000, 3.000000, 0, 1.0};
Point(31) = {2.000000, 3.000000, 0, 1.0};
Point(35) = {3.000000, 3.000000, 0, 1.0};
Line(2) = {1, 2};
Transfinite Line {2} = 1 Using Progression 1;
Line(3) = {2, 3};
Transfinite Line {3} = 1 Using Progression 1;
Line(4) = {3, 4};
Transfinite Line {4} = 1 Using Progression 1;
Line(5) = {4, 1};
Transfinite Line {5} = 1 Using Progression 1;
Line(10) = {2, 6};
Transfinite Line {10} = 1 Using Progression 1;
Line(11) = {6, 7};
Transfinite Line {11} = 1 Using Progression 1;
Line(12) = {7, 3};
Transfinite Line {12} = 1 Using Progression 1;
Line(18) = {6, 10};
Transfinite Line {18} = 1 Using Progression 1;
Line(19) = {10, 11};
Transfinite Line {19} = 1 Using Progression 1;
Line(20) = {11, 7};
Transfinite Line {20} = 1 Using Progression 1;
Line(27) = {3, 15};
Transfinite Line {27} = 1 Using Progression 1;
Line(28) = {15, 16};
Transfinite Line {28} = 1 Using Progression 1;
Line(29) = {16, 4};
Transfinite Line {29} = 1 Using Progression 1;
Line(35) = {7, 19};
Transfinite Line {35} = 1 Using Progression 1;
Line(36) = {19, 15};
Transfinite Line {36} = 1 Using Progression 1;
Line(43) = {11, 23};
Transfinite Line {43} = 1 Using Progression 1;
Line(44) = {23, 19};
Transfinite Line {44} = 1 Using Progression 1;
Line(51) = {15, 27};
Transfinite Line {51} = 1 Using Progression 1;
Line(52) = {27, 28};
Transfinite Line {52} = 1 Using Progression 1;
Line(53) = {28, 16};
Transfinite Line {53} = 1 Using Progression 1;
Line(59) = {19, 31};
Transfinite Line {59} = 1 Using Progression 1;
Line(60) = {31, 27};
Transfinite Line {60} = 1 Using Progression 1;
Line(67) = {23, 35};
Transfinite Line {67} = 1 Using Progression 1;
Line(68) = {35, 31};
Transfinite Line {68} = 1 Using Progression 1;
Line Loop(6) = {2, 3, 4, 5};
Line Loop(14) = {10, 11, 12, -3};
Line Loop(22) = {18, 19, 20, -11};
Line Loop(30) = {-4, 27, 28, 29};
Line Loop(38) = {-12, 35, 36, -27};
Line Loop(46) = {-20, 43, 44, -35};
Line Loop(54) = {-28, 51, 52, 53};
Line Loop(62) = {-36, 59, 60, -51};
Line Loop(70) = {-44, 67, 68, -59};
Plane Surface(7) = {6};
Transfinite Surface {7};
Plane Surface(15) = {14};
Transfinite Surface {15};
Plane Surface(23) = {22};
Transfinite Surface {23};
Plane Surface(31) = {30};
Transfinite Surface {31};
Plane Surface(39) = {38};
Transfinite Surface {39};
Plane Surface(47) = {46};
Transfinite Surface {47};
Plane Surface(55) = {54};
Transfinite Surface {55};
Plane Surface(63) = {62};
Transfinite Surface {63};
Plane Surface(71) = {70};
Transfinite Surface {71};
Physical Surface(30) = {7, 15, 23, 31, 39, 47, 55, 63, 71 };
Physical Surface(31) = {};
Recombine Surface {23};
Recombine Surface {47};
Recombine Surface {71};
Physical Line(72) = {2, 10, 18 };
Physical Line(73) = {19, 43, 67 };
Physical Line(74) = {52, 60, 68 };
Physical Line(75) = {5, 53 };
Physical Line(76) = {29};