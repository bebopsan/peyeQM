// Gmsh project created on Fri Nov 18 10:18:58 2011
Point(1) = {0.57735, 0, 0, 1.0};
Point(2) = {0.28867, 0.5, 0, 1.0};
Point(3) = {-0.28867, 0.5, 0, 1.0};
Point(4) = {-0.57735, 0, 0, 1.0};
Point(5) = {-0.28867, -0.5, 0, 1.0};
Point(6) = {0.28867, -0.5, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line Loop(7) = {2, 3, 4, 5, 6, 1};
Physical Line(8) = {1, 2, 3, 4, 5, 6};

Line(9) = {1, 4};
Line Loop(10) = {9, -3, -2, -1};
Plane Surface(11) = {10};
Line Loop(12) = {9, 4, 5, 6};
Plane Surface(13) = {12};
Physical Surface(14) = {11,13};
