cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {4, 0, 0, cl1};
Point(3) = {4, 2, 0, cl1};
Point(4) = {0, 2, 0, cl1};
Line(1) = {1, 2};
Transfinite Line {1} = 3 Using Progression 1;
Line(2) = {2, 3};
Transfinite Line {2} = 2 Using Progression 1;
Line(3) = {3, 4};
Transfinite Line {3} = 3 Using Progression 1;
Line(4) = {4, 1};
Transfinite Line {4} = 2 Using Progression 1;
Line Loop(6) = {3, 4, 1, 2};
Plane Surface(6) = {6};
Transfinite Surface {6} Alternated;
Recombine Surface {6};
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};
Physical Surface(11) = {6};