Point(1) = {0.000000, 0.000000, 0, 1.0};
Point(2) = {2.000000, 0.000000, 0, 1.0};
Point(3) = {2.000000, 2.000000, 0, 1.0};
Point(4) = {0.000000, 2.000000, 0, 1.0};
Point(6) = {4.000000, 0.000000, 0, 1.0};
Point(7) = {4.000000, 2.000000, 0, 1.0};
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
Line Loop(6) = {2, 3, 4, 5};
Line Loop(14) = {10, 11, 12, -3};
Plane Surface(7) = {6};
Transfinite Surface {7};
Plane Surface(15) = {14};
Transfinite Surface {15};
Physical Surface(30) = {7, 15 };
Recombine Surface {15};
Recombine Surface {7};
Physical Line(16) = {2, 10 };
Physical Line(17) = {11 };
Physical Line(18) = {4, 12 };
Physical Line(19) = {5 };
