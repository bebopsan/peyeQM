Point(1) = {0.000000, 0.000000, 0, 1.0};
Point(2) = {2.000000, 0.000000, 0, 1.0};
Point(3) = {2.000000, 2.000000, 0, 1.0};
Point(4) = {0.000000, 2.000000, 0, 1.0};
Point(6) = {4.000000, 0.000000, 0, 1.0};
Point(7) = {4.000000, 2.000000, 0, 1.0};
Point(10) = {6.000000, 0.000000, 0, 1.0};
Point(11) = {6.000000, 2.000000, 0, 1.0};
Point(14) = {8.000000, 0.000000, 0, 1.0};
Point(15) = {8.000000, 2.000000, 0, 1.0};
Point(18) = {10.000000, 0.000000, 0, 1.0};
Point(19) = {10.000000, 2.000000, 0, 1.0};
Point(22) = {12.000000, 0.000000, 0, 1.0};
Point(23) = {12.000000, 2.000000, 0, 1.0};
Point(26) = {14.000000, 0.000000, 0, 1.0};
Point(27) = {14.000000, 2.000000, 0, 1.0};
Point(30) = {16.000000, 0.000000, 0, 1.0};
Point(31) = {16.000000, 2.000000, 0, 1.0};
Point(34) = {18.000000, 0.000000, 0, 1.0};
Point(35) = {18.000000, 2.000000, 0, 1.0};
Line(2) = {1, 2};
Transfinite Line {2} = 2 Using Progression 1;
Line(3) = {2, 3};
Transfinite Line {3} = 2 Using Progression 1;
Line(4) = {3, 4};
Transfinite Line {4} = 2 Using Progression 1;
Line(5) = {4, 1};
Transfinite Line {5} = 2 Using Progression 1;
Line(10) = {2, 6};
Transfinite Line {10} = 2 Using Progression 1;
Line(11) = {6, 7};
Transfinite Line {11} = 2 Using Progression 1;
Line(12) = {7, 3};
Transfinite Line {12} = 2 Using Progression 1;
Line(18) = {6, 10};
Transfinite Line {18} = 2 Using Progression 1;
Line(19) = {10, 11};
Transfinite Line {19} = 2 Using Progression 1;
Line(20) = {11, 7};
Transfinite Line {20} = 2 Using Progression 1;
Line(26) = {10, 14};
Transfinite Line {26} = 2 Using Progression 1;
Line(27) = {14, 15};
Transfinite Line {27} = 2 Using Progression 1;
Line(28) = {15, 11};
Transfinite Line {28} = 2 Using Progression 1;
Line(34) = {14, 18};
Transfinite Line {34} = 2 Using Progression 1;
Line(35) = {18, 19};
Transfinite Line {35} = 2 Using Progression 1;
Line(36) = {19, 15};
Transfinite Line {36} = 2 Using Progression 1;
Line(42) = {18, 22};
Transfinite Line {42} = 2 Using Progression 1;
Line(43) = {22, 23};
Transfinite Line {43} = 2 Using Progression 1;
Line(44) = {23, 19};
Transfinite Line {44} = 2 Using Progression 1;
Line(50) = {22, 26};
Transfinite Line {50} = 2 Using Progression 1;
Line(51) = {26, 27};
Transfinite Line {51} = 2 Using Progression 1;
Line(52) = {27, 23};
Transfinite Line {52} = 2 Using Progression 1;
Line(58) = {26, 30};
Transfinite Line {58} = 2 Using Progression 1;
Line(59) = {30, 31};
Transfinite Line {59} = 2 Using Progression 1;
Line(60) = {31, 27};
Transfinite Line {60} = 2 Using Progression 1;
Line(66) = {30, 34};
Transfinite Line {66} = 2 Using Progression 1;
Line(67) = {34, 35};
Transfinite Line {67} = 2 Using Progression 1;
Line(68) = {35, 31};
Transfinite Line {68} = 2 Using Progression 1;
Line Loop(6) = {2, 3, 4, 5};
Line Loop(14) = {10, 11, 12, -3};
Line Loop(22) = {18, 19, 20, -11};
Line Loop(30) = {26, 27, 28, -19};
Line Loop(38) = {34, 35, 36, -27};
Line Loop(46) = {42, 43, 44, -35};
Line Loop(54) = {50, 51, 52, -43};
Line Loop(62) = {58, 59, 60, -51};
Line Loop(70) = {66, 67, 68, -59};
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
Physical Surface(100) = {7, 15, 23, 31, 39, 47, 55, 63, 71 };
Recombine Surface {71};
Physical Line(72) = {2, 10, 18, 26, 34, 42, 50, 58, 66 };
Physical Line(73) = {67 };
Physical Line(74) = {4, 12, 20, 28, 36, 44, 52, 60, 68 };
Physical Line(75) = {5 };
