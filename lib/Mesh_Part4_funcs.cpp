#include "Meshfuncs.h"
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <cfloat>
#include <time.h>
#include <engine.h>

#pragma comment ( lib, "libmat.lib" )
#pragma comment ( lib, "libmx.lib" )
#pragma comment ( lib, "libmex.lib" )
#pragma comment ( lib, "libeng.lib" )

using namespace std;

void nurb::create_xmsh(string filename, int degree) 
{
	ofstream xmsh_file;
	xmsh_file.open(filename + "_final.xmsh");
	
	// get system time
	time_t cur_time = time(NULL);
	char str[26];
	ctime_s(str, sizeof str, &cur_time);


	xmsh_file << "% MATLAB output file for " << filename << endl;
	xmsh_file << "% Completed at " << str << endl;
	xmsh_file << "% Degree of mesh:" << endl << degree << endl << endl;

	// loop through the node list

	xmsh_file << "% Node List" << endl << "% NUMBER OF NODES" << endl << node_list.size() << endl;
	unsigned int node_size = unsigned int(node_list.size());
	for (unsigned int i = 0; i < node_size; i++) {
		xmsh_file << i+1 << "," << node_list[i][0] << "," << node_list[i][1] << "," << node_list[i][2] << endl;
	}
	xmsh_file << "% END OF NODE LIST " << endl << endl;

	// loop through the connectivity array

	xmsh_file << "% Connectivity array" << endl << "% NUMBER OF TRIANGLES, NODES PER TRIANGLE" << endl << triangles.size() << "," << nodes_in_triangle << endl;
	unsigned int tri_size = unsigned (triangles.size());
	for (unsigned int i = 0; i < tri_size; i++) {
		xmsh_file << i + 1 << ",";
		for (int j = 0; j < nodes_in_triangle; j++) {
			if (j == nodes_in_triangle - 1)
				xmsh_file << triangles[i]->controlP[j] + 1 << endl;
			else
				xmsh_file << triangles[i]->controlP[j] + 1 << ",";
		}
	}
	xmsh_file << "% END OF CONNECTIVITY ARRAY" << endl << endl;

	// now give the NURBS curve data, THIS IS TEMPORARY
	xmsh_file << " % NURBS Curve Information" << endl;
	xmsh_file << " % The format is: " << endl << " % NURBS CURVE,ELEM NUM,XPOINT,YPOINT,WEIGHT" << endl;
	for (int i = 0; i < num_curves; i++) {
		for (int j = 0; j < Elem_list[i]->n_el; j++) {
			for (int k = 0; k < degree; k++) {
				xmsh_file << i << "," << j << "," << Elem_list[i]->elem_Geom[j].controlP[k][0] << "," << Elem_list[i]->elem_Geom[j].controlP[k][1] << "," << Elem_list[i]->elem_Geom[j].controlP[k][2] << endl;
			}
		}
	}
	xmsh_file << " % END OF CURVE INFORMATION" << endl << endl;

	// now give the data that connects the nurbs curves to the triangles
	xmsh_file << " % Information connecting NURBS curve and Triangles" << endl << tri_NURB_elem_section_side.size() << endl;

	for (unsigned int i = 0; i < int(tri_NURB_elem_section_side.size()); i++) {
		vector<unsigned int> tri_row = tri_NURB_elem_section_side[i];
		xmsh_file << i + 1 << "," << tri_row[0] << "," << tri_row[1] << "," << tri_row[2] << "," << tri_row[4] << "," << tri_row[5] << endl;
	}
	xmsh_file << " % END OF CONNECTING INFORMATION" << endl;
	xmsh_file.close();
}

void nurb::display_mesh(string filename)
{
	Engine *mat_eng;
	mat_eng = engOpen("null");
	engSetVisible(mat_eng, false);
	engEvalString(mat_eng, "close all");
	//engEvalString(mat_eng, "figure");

	// for now I am just going to display the control points
	// first I need to convert the node_list from a 2D vector into a 2D array

	const int num_nodes = int(node_list.size());
	double **node_array = vec_to_array(node_list, num_nodes);
	double *column1 = node_array[0];
	double *column2 = node_array[1];

	mxArray *mat_nodes1 = mxCreateDoubleMatrix(num_nodes, 1, mxREAL);
	mxArray *mat_nodes2 = mxCreateDoubleMatrix(num_nodes, 1, mxREAL);
	// copy over the data
	memcpy((void *)mxGetPr(mat_nodes1), (void *)column1, sizeof(double)*num_nodes);
	memcpy((void *)mxGetPr(mat_nodes2), (void *)column2, sizeof(double)*num_nodes);
	engPutVariable(mat_eng, "x_coor", mat_nodes1);
	engPutVariable(mat_eng, "y_coor", mat_nodes2);


	engEvalString(mat_eng, "scatter(x_coor, y_coor, '.'); hold on");
	engEvalString(mat_eng, "title('Mesh', 'Interpreter','latex','FontSize',16)");
	engEvalString(mat_eng, "axis equal");

	mxDestroyArray(mat_nodes1);
	mxDestroyArray(mat_nodes2);

	// now loop through all of the global edges and display the connecting nurbs curves
	// first I need to create a template for the basis evaluations since this can be reused from element to element
	vector<vector<double>> N(11, vector<double>(degree + 1, 0.0));
	for (int j = 0; j < 11; j++) {
		for (int i = 0; i <= degree; i++) {
			double xi = double(j) / 10.0;
			N[j][i] = n_choose_k(degree, i) * pow(xi, i) * pow(1 - xi, degree - i);
		}
	}
	mxArray *x_points = mxCreateDoubleMatrix(11, 1, mxREAL);  // create the mxArray of the correct size
	mxArray *y_points = mxCreateDoubleMatrix(11, 1, mxREAL);
	for (size_t i = 0; i < global_edges.size(); i++) {
		vector<vector<double>> R = eval_edges(N, int(i));
			for (int j = 0; j < 11; j++) { // loop through all of the xi values
				column1[j] = 0;
				column2[j] = 0;
				for (int k = 0; k <= degree; k++) {
					column1[j] += node_list[global_edges[i][k]][0] * R[j][k];
					column2[j] += node_list[global_edges[i][k]][1] * R[j][k];
				}
			}

		memcpy((void *)mxGetPr(x_points), (void *)column1, sizeof(double) * 11);   // copy the data over
		memcpy((void *)mxGetPr(y_points), (void *)column2, sizeof(double) * 11);
		engPutVariable(mat_eng, "x_coor", x_points);   // set the variable in the matlab workspace
		engPutVariable(mat_eng, "y_coor", y_points);
		engEvalString(mat_eng, "plot(x_coor, y_coor, 'r')");
		engEvalString(mat_eng, "grid on; grid minor");

	}
	// free all the memory
	mxDestroyArray(x_points);
	mxDestroyArray(y_points);
	delete[] column1;
	delete[] column2;
	delete[] node_array;
	//engClose(mat_eng);
}

double** nurb::vec_to_array(vector<vector<double>> &vec, unsigned int rows)
{
	double** temp;
	temp = new double*[2];
	temp[0] = new double [rows];
	temp[1] = new double [rows];
	for (unsigned int i = 0; i < rows; i++) {
		temp[0][i] = vec[i][0];
		temp[1][i] = vec[i][1];

	}
	return temp;
}


vector<vector<double>> nurb::eval_edges(vector<vector<double>> N, int edge)
{
	// First I need to figure out the weight to determine the rational basis function evaluations
	vector<vector<double>> R(11, vector<double>(degree + 1, 0.0));
	for (int i = 0; i < 11; i++) {
		double w_tot = 0;
		for (int j = 0; j <= degree; j++) {
			w_tot += node_list[global_edges[edge][j]][2] * N[i][j];
		}
		for (int j = 0; j <= degree; j++) {
			R[i][j] = node_list[global_edges[edge][j]][2] * N[i][j] / w_tot;
			if (isnan(R[i][j]))
				R[i][j] = 0.0;
		}
	}
	return R;
}