// Gmsh project created on Sat Apr  6 20:16:57 2013
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0.92387953251, 0.38268343236, 0, 1.0};
Point(4) = {5, 0, 0, 1.0};
Point(5) = {4.61939766256, 1.91341716183, 0, 1.0};
Circle(1) = {3, 1, 2};
Line(2) = {2, 4};
Circle(3) = {4, 1, 5};
Line(4) = {5, 3};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Transfinite Line {1} = 10 Using bump 1;
Transfinite Line {2, -4} = 50 Using Progression 1.08;
Transfinite Line {3} = 10 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};
Physical Line(7) = {1};
Physical Line(8) = {3};
Physical Line(9) = {2, -4};
Physical Surface(10) = {6};
