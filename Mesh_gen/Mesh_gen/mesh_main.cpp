#include <iostream>
#include <stdlib.h>
#include "Meshfuncs.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


int main(int argc, char*argv[]) {
	nurb *action = new nurb;
	action->fact_table();
	cout << "Beginning Initialization of Mesh.........." << endl;

	//  Need to do a readin file parser for the .spline file
	// the first command line argument should be the filename
	action->head_nurb = action->readSpline(argv[1]);
	action->NURBS2poly(action->head_nurb);
	// change eval_bern to change the evaluation rules.  right now it does it at xi = 0, 0.5, 1

	// Part 2 Making the Linear Mesh
	action->create_geo_file(argv[1]);
	action->call_gmsh(argv[1]);
	action->readMsh(argv[1]);

	// Part 3 Smoothing and Degree Elevating the mesh
	action->smoothMesh();

	
	return 0;
}