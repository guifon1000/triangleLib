// Inputs
squareSide = 200; //m
meshThickness = squareSide / 10; 
N1=5;
N2=4;

// Geometry
Point(1) = {-0.25,2.,0.15};
Point(2) = {0.3,2.,-0.45};
Point(3) = {0.5,-2.,-0.3};
Point(4) = {-1.,-2.,0.35};
Line(1) = {1, 2};				// bottom line
Line(2) = {2, 3};				// right line
Line(3) = {3, 4};				// top line
Line(4) = {4, 1};				// left line
Line Loop(5) = {1, 2, 3, 4}; 	
Ruled Surface(6) = {5};

Transfinite Line {4, 2} = N1 Using Progression 1;

Transfinite Line {1, 3} = N2 Using Progression 1;
//Transfinite surface:
Transfinite Surface {6};
Recombine Surface {6};
 
