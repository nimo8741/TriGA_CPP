//  This is for NURBS Curve 1
Point(1) = {1.000000, 0.000000, 0.0, 1.0};
Point(2) = {0.929788, 0.368095, 0.0, 1.0};
Point(3) = {0.707107, 0.707107, 0.0, 1.0};
Point(4) = {0.368095, 0.929788, 0.0, 1.0};
Point(5) = {0.000000, 1.000000, 0.0, 1.0};
Point(6) = {-0.368095, 0.929788, 0.0, 1.0};
Point(7) = {-0.707107, 0.707107, 0.0, 1.0};
Point(8) = {-0.929788, 0.368095, 0.0, 1.0};
Point(9) = {-1.000000, 0.000000, 0.0, 1.0};
Point(10) = {-0.929788, -0.368095, 0.0, 1.0};
Point(11) = {-0.707107, -0.707107, 0.0, 1.0};
Point(12) = {-0.368095, -0.929788, 0.0, 1.0};
Point(13) = {0.000000, -1.000000, 0.0, 1.0};
Point(14) = {0.368095, -0.929788, 0.0, 1.0};
Point(15) = {0.707107, -0.707107, 0.0, 1.0};
Point(16) = {0.929788, -0.368095, 0.0, 1.0};

Line(1) = {1, 2};
Physical Line("Knot Vector Section 1") = {1, 2};
Line(2) = {2, 3};
Physical Line("Knot Vector Section 2") = {2, 3};
Line(3) = {3, 4};
Physical Line("Knot Vector Section 3") = {3, 4};
Line(4) = {4, 5};
Physical Line("Knot Vector Section 4") = {4, 5};
Line(5) = {5, 6};
Physical Line("Knot Vector Section 5") = {5, 6};
Line(6) = {6, 7};
Physical Line("Knot Vector Section 6") = {6, 7};
Line(7) = {7, 8};
Physical Line("Knot Vector Section 7") = {7, 8};
Line(8) = {8, 9};
Physical Line("Knot Vector Section 8") = {8, 9};
Line(9) = {9, 10};
Physical Line("Knot Vector Section 9") = {9, 10};
Line(10) = {10, 11};
Physical Line("Knot Vector Section 10") = {10, 11};
Line(11) = {11, 12};
Physical Line("Knot Vector Section 11") = {11, 12};
Line(12) = {12, 13};
Physical Line("Knot Vector Section 12") = {12, 13};
Line(13) = {13, 14};
Physical Line("Knot Vector Section 13") = {13, 14};
Line(14) = {14, 15};
Physical Line("Knot Vector Section 14") = {14, 15};
Line(15) = {15, 16};
Physical Line("Knot Vector Section 15") = {15, 16};
Line(16) = {16, 1};
Physical Line("Knot Vector Section 16") = {16, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};


//  This is for NURBS Curve 2
Point(17) = {4.000000, 4.000000, 0.0, 1.0};
Point(18) = {-4.000000, 4.000000, 0.0, 1.0};
Point(19) = {-4.000000, -4.000000, 0.0, 1.0};
Point(20) = {4.000000, -4.000000, 0.0, 1.0};

Line(17) = {17, 18};
Physical Line("Knot Vector Section 17") = {17, 18};
Line(18) = {18, 19};
Physical Line("Knot Vector Section 18") = {18, 19};
Line(19) = {19, 20};
Physical Line("Knot Vector Section 19") = {19, 20};
Line(20) = {20, 17};
Physical Line("Knot Vector Section 20") = {20, 17};

Line Loop(2) = {17, 18, 19, 20};


Plane Surface(1) = {1, 2};
Physical Surface("Domain of interest") = {1};
