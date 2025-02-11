// Gmsh project created on Fri May 17 18:23:14 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {0, 0, 1, 1.0};
//+
Point(6) = {0, 1, 1, 1.0};
//+
Point(7) = {1, 1, 1, 1.0};
//+
Point(8) = {1, 0, 1, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {5, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 6};
//+
Line(8) = {6, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {4, 8};
//+
Line(11) = {3, 7};
//+
Line(12) = {2, 6};
//+
Circle(13) = {0.5, 0, 0.5, 0.4, 0, 2*Pi};
//+
Circle(14) = {0.5, 1, 0.5, 0.4, 0, 2*Pi};
//+
Circle(15) = {0.5, 0.5, 0, 0.2, 0, 2*Pi};
//+
Circle(16) = {0.5, 0.5, 1, 0.2, 0, 2*Pi};
//+
Rotate {{1, 0, 0}, {0.5, 0, 0.5}, Pi/2} {
  Curve{13}; 
}
//+
Rotate {{1, 0, 0}, {0.5, 1, 0.5}, Pi/2} {
  Curve{14}; 
}
//+
Curve Loop(1) = {13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16};
//+
Curve Loop(4) = {15};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {16};
//+
Plane Surface(4) = {5};
//+
Curve Loop(6) = {1, 10, -5, -9};
//+
Curve Loop(7) = {13};
//+
Plane Surface(5) = {6, 7};
//+
Curve Loop(8) = {3, 12, -7, -11};
//+
Curve Loop(9) = {14};
//+
Plane Surface(6) = {8, 9};
//+
Curve Loop(10) = {1, 2, 3, 4};
//+
Curve Loop(11) = {15};
//+
Plane Surface(7) = {10, 11};
//+
Curve Loop(12) = {5, 6, 7, 8};
//+
Curve Loop(13) = {16};
//+
Plane Surface(8) = {12, 13};
//+
Curve Loop(14) = {10, 6, -11, -2};
//+
Plane Surface(9) = {14};
//+
Curve Loop(15) = {4, 9, -8, -12};
//+
Plane Surface(10) = {15};
//+
Surface Loop(1) = {1, 5, 7, 9, 8, 6, 10, 2, 4, 3};
//+
Volume(1) = {1};
