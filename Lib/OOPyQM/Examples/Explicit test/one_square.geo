Point(1) = {0.000000, 0.000000, 0, 1.0};
Point(2) = {2.000000, 0.000000, 0, 1.0};
Point(3) = {2.000000, 2.000000, 0, 1.0};
Point(4) = {0.000000, 2.000000, 0, 1.0};
Line(2) = {1, 2};
Transfinite Line {2} = 2 Using Progression 1;
Line(3) = {2, 3};
Transfinite Line {3} = 2 Using Progression 1;
Line(4) = {3, 4};
Transfinite Line {4} = 2 Using Progression 1;
Line(5) = {4, 1};
Transfinite Line {5} = 2 Using Progression 1;
Line Loop(6) = {2, 3, 4, 5};
Plane Surface(7) = {6};
Transfinite Surface {7};
Physical Surface(30) = {7 };
Recombine Surface {7};
Physical Line(8) = {2 };
Physical Line(9) = {3 };
Physical Line(10) = {4 };
Physical Line(11) = {5 };
