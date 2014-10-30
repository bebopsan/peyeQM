// Gmsh project created on Tue Apr  9 11:45:19 2013
Point(1) = {2, 0, 0, 1.0};
Point(2) = {3, 1, 0, 1.0};
Point(3) = {3, 2, 0, 1.0};
Point(4) = {2, 2, 0, 1.0};
Point(5) = {1, 2, 0, 1.0};
Point(6) = {1, 3, 0, 1.0};
Point(7) = {0, 2, 0, 1.0};
Point(8) = {0.5, 1.5, 0, 1.0};
Point(9) = {1, 1, 0, 1.0};
Point(10) = {2, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};
Line(11) = {10, 3};
Line(12) = {9, 4};
Line(13) = {8, 5};
Line Loop(14) = {10, 1, 2, -11};
Plane Surface(15) = {14};
Line Loop(16) = {9, 11, 3, -12};
Plane Surface(17) = {16};
Line Loop(18) = {8, 12, 4, -13};
Plane Surface(19) = {18};
Line Loop(20) = {7, 13, 5, 6};
Plane Surface(21) = {20};
Physical Line(22) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Physical Line(23) = {11, 12, 13};
Physical Surface(24) = {15, 17, 19, 21};
Transfinite Line {1, 2, 10, 11, 9, 12, 3, 4, 8, 13, 5, 6, 7} = 5 Using Progression 1;
Transfinite Surface {15};
Transfinite Surface {17};
Transfinite Surface {21};
Transfinite Surface {19};
Recombine Surface {15, 17, 19, 21};