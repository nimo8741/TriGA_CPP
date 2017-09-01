#include "Meshfuncs.h"
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <time.h>

using namespace std;

void nurb::create_xmsh(string filename, int degree)
{
	ofstream xmsh_file;
	xmsh_file.open(path_to_file + "\\IO_files\\xmsh_files\\" + filename + ".xmsh");

	// get system time
	time_t cur_time = time(NULL);
	char str[26];
	#ifdef _WIN32
		ctime_s(str, sizeof str, &cur_time);
	#else
		sprintf(str, "%s", ctime(&cur_time));
	#endif


	xmsh_file << "% XMSH output file for " << filename << endl;
	xmsh_file << "% Completed at " << str << endl;
	xmsh_file << "% Degree of mesh:" << endl << degree << endl << endl;

	// loop through the node list

	xmsh_file << "% Node List" << endl << "% NUMBER OF NODES" << endl << node_list.size() << endl;
	unsigned int node_size = int(node_list.size());
	for (unsigned int i = 0; i < node_size; i++) {
		xmsh_file << i + 1 << "," << node_list[i][0] << "," << node_list[i][1] << "," << node_list[i][2] << endl;
	}
	xmsh_file << "% END OF NODE LIST " << endl << endl;

	// loop through the connectivity array

	xmsh_file << "% Connectivity array" << endl << "% NUMBER OF TRIANGLES, NODES PER TRIANGLE" << endl << triangles.size() << "," << nodes_in_triangle << endl;
	unsigned int tri_size = unsigned(triangles.size());
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

	// loop through global edges list
	unsigned int edge_size = unsigned(global_edges.size());
	xmsh_file << "% Global Edge List" << "% NUMBER OF EDGES, NODES PER EDGE" << endl << edge_size << degree + 1 << endl;
	for (unsigned int i = 0; i < edge_size; i++) {
		xmsh_file << i << ",";
		for (int j = 0; j < degree; j++) {
			xmsh_file << global_edges[i][j] << ",";
		}
		xmsh_file << global_edges[i][degree] << endl;
	}
	xmsh_file << "% END OF GLOBAL EDGES" << endl << endl;
	xmsh_file.close();
}

void nurb::display_mesh(string filename)
{
	ofstream control_file;
	control_file.open(path_to_file + "\\IO_files\\dat_files\\control.dat");

	double xmin = node_list[0][0];
	double xmax = node_list[0][0];
	double ymin = node_list[0][1];
	double ymax = node_list[0][1];

	control_file << "# This is control.dat" << endl;
	for (unsigned int i = 0; i < node_list.size(); i++) {
		control_file << node_list[i][0] << " " << node_list[i][1] << endl;

		// update the min and max for sizing the plot
		if (node_list[i][0] < xmin)
			xmin = node_list[i][0];
		else if (node_list[i][0] > xmax)
			xmax = node_list[i][0];
		if (node_list[i][1] < ymin)
			ymin = node_list[i][1];
		else if (node_list[i][1] > ymax)
			ymax = node_list[i][1];
	}
	control_file.close();
	int xmin_i = (xmin < 0) ? int(floor(xmin)) : int(ceil(xmin));
	int xmax_i = (xmax < 0) ? int(floor(xmax)) : int(ceil(xmax));
	int ymin_i = (ymin < 0) ? int(floor(ymin)) : int(ceil(ymin));
	int ymax_i = (ymax < 0) ? int(floor(ymax)) : int(ceil(ymax));

	vector<vector<double>> N(11, vector<double>(degree + 1, 0.0));
	for (int j = 0; j < 11; j++) {
		for (int i = 0; i <= degree; i++) {
			double xi = double(j) / 10.0;
			N[j][i] = n_choose_k(degree, i) * pow(xi, i) * pow(1 - xi, degree - i);
		}
	}

	ofstream edge_file;
	edge_file.open(path_to_file + "\\IO_files\\dat_files\\edges.dat");
	edge_file << "# This is edges.dat" << endl;
	for (size_t i = 0; i < global_edges.size(); i++) {
		vector<vector<double>> R = eval_edges(N, int(i));
		vector<vector<double>> points(11, vector<double>(2, 0));
		for (int j = 0; j < 11; j++) { // loop through all of the xi values
			for (int k = 0; k <= degree; k++) {
				points[j][0] += node_list[global_edges[i][k]][0] * R[j][k];
				points[j][1] += node_list[global_edges[i][k]][1] * R[j][k];
			}
			edge_file << points[j][0] << " " << points[j][1] << endl;
		}
		edge_file << endl;

	}
	string cmd = "\"" + path_to_file + "gnuplot\\build\\gnuplot.exe\"";
	string args = "\"" + path_to_file + "plot_cmds.plot\"";

	// now I need to create a verison of path_to_file which contains double backslashes instead of single ones
	string d_back_path = path_to_file;
	for (unsigned int i = 0; i < d_back_path.size(); i++) {
		if (d_back_path[i] == '\\') { // there is a single backslash
			d_back_path.insert(i, "\\");
			i++;
		}
	}
	// now I need to rewrite the plot_cmds.plot
	ofstream plot_file;
	plot_file.open(path_to_file + "\\plot_cmds.plot");
	plot_file << "cd '" + d_back_path + "IO_files\\\\dat_files'" << endl;
	plot_file << "set size square" << endl << "plot 'control.dat' with points linecolor rgb 'blue' pointtype 7 ps 0.2 notitle, \\" << endl;
	plot_file << "     'edges.dat' with lines linecolor rgb 'red' lw 0.5 notitle" << endl << "set size square" << endl << "pause -1 \"Hit return to continue\"";
	plot_file.close();


	const int size = 1000;
	char runline[size];
	
	#ifdef _WIN32
		sprintf_s(runline, size, "\"%s %s\"", cmd.c_str(), args.c_str());
	#else
		sprintf(runline, "\"%s %s\"", cmd.c_str(), args.c_str());
	#endif
	
	system((char *)runline);   // this line runs gmsh with default settings

}

double** nurb::vec_to_array(vector<vector<double>> &vec, unsigned int rows)
{
	double** temp;
	temp = new double*[2];
	temp[0] = new double[rows];
	temp[1] = new double[rows];
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
