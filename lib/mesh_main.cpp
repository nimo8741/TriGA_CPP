#include <iostream>
#include <stdlib.h>
#include "Meshfuncs.h"
#include <Eigen/Dense>
#include <cstdarg>

using namespace std;
using namespace Eigen;


int main(int argc, char*argv[]) {
	if (argc == 1) {
		cout << "Please provide input arguments" << endl;
		return 0;
	}
	nurb *action = new nurb;

	int degree;
	if (argc == 3) {
		degree = stoi(argv[2]);
	}
	else
		degree = 3;   //default to order 3

	action->fact_table(degree);
	//  Need to do a readin file parser for the .spline file
	// the first command line argument should be the filename
	action->head_nurb = action->readSpline(argv[1]);
	action->NURBS2poly(action->head_nurb);

	// Part 2 Making the Linear Mesh
	action->create_geo_file(argv[1]);
	action->call_gmsh(argv[1]);


	action->readMsh(argv[1], degree);
	action->create_side_nodes();

	// Part 3 Smoothing and Degree Elevating the mesh
	action->smoothMesh(degree);

	// Part 4 Write output file to matlab to then be diplayed
	action->create_xmsh(argv[1], degree);
	action->display_mesh(argv[1]);


	return 0;
}
