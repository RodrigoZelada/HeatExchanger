// Gmsh project created on Mon Aug  5 18:19:37 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 0, 1, 1.0};
//+
Point(6) = {1, 0, 1, 1.0};
//+
Point(7) = {1, 1, 1, 1.0};
//+
Point(8) = {0, 1, 1, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Circle(13) = {0.5, 0, 0.5, 0.1, 0, 2*Pi};
//+
Circle(14) = {0.5, 1, 0.5, 0.1, 0, 2*Pi};
//+
Circle(15) = {0.5, 0.25, 0, 0.1, 0, 2*Pi};
//+
Circle(16) = {0.5, 0.75, 1, 0.1, 0, 2*Pi};
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
Curve Loop(3) = {15};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {1, 10, -5, -9};
//+
Curve Loop(6) = {13};
//+
Plane Surface(5) = {5, 6};
//+
Curve Loop(7) = {1, 2, 3, 4};
//+
Curve Loop(8) = {15};
//+
Plane Surface(6) = {7, 8};
//+
Curve Loop(9) = {3, 12, -7, -11};
//+
Curve Loop(10) = {14};
//+
Plane Surface(7) = {9, 10};
//+
Curve Loop(11) = {6, 7, 8, 5};
//+
Curve Loop(12) = {16};
//+
Plane Surface(8) = {11, 12};
//+
Curve Loop(13) = {10, 6, -11, -2};
//+
Plane Surface(9) = {13};
//+
Curve Loop(14) = {9, -8, -12, 4};
//+
Plane Surface(12) = {14};
//+
Surface Loop(1) = {5, 6, 9, 8, 7, 12, 2, 4, 3, 1};
//+
Volume(1) = {1};
