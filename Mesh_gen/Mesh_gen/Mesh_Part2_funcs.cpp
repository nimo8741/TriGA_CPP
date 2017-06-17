#include "Meshfuncs.h"
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


/**********************************************************************************
Function prototype:
void nurb::create_geo_file(string filename)

Function description:
This function ceates the .geo file which would then be passed into the linear mesh
generator gmsh


Precondition:
This function requires the name of the file to be created.  This is the same filename
used to read the .spline file.

Postcondition:
This function writes the .geo file and has the filename as the same as the original
filename passed in at the beginning of this program's execution
**********************************************************************************/

void nurb::create_geo_file(string filename)
{
	// First and foremost, create the .geo file
	ofstream geo_file;
	geo_file.open(filename + ".geo");

	int p_count = 1;   // point counter
	int l_count = 1;   // line (edge) counter
	const int size = 200;
	char line [size];

	// loop through all of the NURBS curves
	for (unsigned int i = 0; i < Elem_list.size(); i++) {
		int p_start = p_count;
		int l_start = l_count;
		Bezier_handle * current = Elem_list[i];
		geo_file << "//  This is for NURBS Curve " << to_string(i + 1) << endl;

		// first loop through all of the quadrature points
		for (int j = 0; j < current->n_el; j++) {                 // loop through each of the elements
			Bezier_elem cur_elem = current->elem_Geom[j];

			for (unsigned int k = 0; k < current->xi_evals.size(); k++) {        // loop through each of the quadrature points within the element
				if (k != current->xi_evals.size() - 1) {
					sprintf_s(line, size, "Point(%d) = {%f, %f, 0.0, 1.0};", p_count, cur_elem.xi_real_loc[k][0], cur_elem.xi_real_loc[k][1]);
					geo_file << line << endl;
					p_count++;
				}
				
			}
		}
		// now loop through all of the connecting edges
		geo_file << endl;
		for (int j = 0; j < p_count - p_start; j++) {
			if (j == p_count - p_start - 1) {
				sprintf_s(line, size, "Line(%d) = {%d, %d};", l_count, j + p_start, p_start);

			}
			else {
				sprintf_s(line, size, "Line(%d) = {%d, %d};", l_count, j + p_start, j + 1 + p_start);
			}

			geo_file << line << endl;
			l_count++;
		}
		geo_file << endl;
		// now need to add the line loop
		string line_loop = "Line Loop(";
		line_loop += to_string(i + 1);
		line_loop += ") = {";
		for (int j = l_start; j < l_count; j++) {
			line_loop += to_string(j);
			line_loop += ", ";
		}
		line_loop.replace(line_loop.end() - 2, line_loop.end(), "};");
		geo_file << line_loop << endl;

		// now make the entry for the physical line this needs the reference the boundary conditions from the .spline file
		Geom_data * cur_Nurb_geo = head_nurb;
		
		// this loop is so that it goes to the correct NURBS curve data
		for (unsigned int next = 0; next < i; next++) {
			cur_Nurb_geo = cur_Nurb_geo->next;
		}

		vector<vector<int>> boundary_group (1, vector<int>(1,0));   // this initializes a dynamic 2D vector with a single 1x1 entry contains a 0 for element 0

		// each row in boundary_group corresponds to a different boundary condition.  These start at condition 0 and increase by 1
		// each value in it denotes the index of the element which belongs to that boundary condition

		// now I need to loop through all of the other elements and find which boundary condition applies to them by which side the element lives on
		for (int j = 1; j < current->n_el; j++) {
			int side = current->elem_Geom[j].side;
			int b_index = cur_Nurb_geo->B_flag[side];

			// now loop through boundary group to 
			bool found = false;
			for (unsigned int k = 0; k < boundary_group.size(); k++) {
				if (cur_Nurb_geo->B_flag[current->elem_Geom[j].side] == cur_Nurb_geo->B_flag[current->elem_Geom[k].side]) {
					boundary_group[k].push_back(j);
					found = true;
				}
			}
			if (found == false) {
				boundary_group.push_back(vector<int>(1, j));

			}
		}

		// now I need to translate these elements into the lines which were used previously. I know that there will be p number of lines in each element
		for (unsigned int j = 1; j <= boundary_group.size(); j++) {
			string phy_line = "Physical Line(\"Boundary Group ";
			phy_line += to_string(i + j);
			phy_line += "\") = {";
			for (unsigned int k = 0; k < boundary_group[j - 1].size(); k++) {
				phy_line += to_string(((boundary_group[j - 1][k] + 1) << 1) + l_start - 2);
				phy_line += ", ";
				phy_line += to_string(((boundary_group[j - 1][k] + 1) << 1) + l_start - 1);
				phy_line += ", ";

			}
			phy_line.replace(phy_line.end() - 2, phy_line.end(), "};");
			geo_file << phy_line << endl;

		}


		geo_file << endl << endl;
	}

	// now add all of these loops into a plane surface
	string plane = "Plane Surface(1) = {";    // there is only going to be a single plane surface
	for (unsigned int i = 1; i <= Elem_list.size(); i++) {
		plane += to_string(i);
		plane += ", ";
	}
	plane.replace(plane.end() - 2, plane.end(), "};");
	geo_file << plane << endl;

	geo_file << "Physical Surface(\"Domain of interest\") = {1};" << endl;

	geo_file.close();
}


/**********************************************************************************
Function prototype:
void nurb::call_gmsh(string filename

Function description:
This function calls the linear mesher gmsh using the previously created .geo file

Precondition:
This function uses the filename used in the creation of the .geo file

Postcondition:
After initiating the linear mesh generator, execution of the program is passed off
for the time being and the program awaits for gmsh to create its .msh output file
**********************************************************************************/

void nurb::call_gmsh(string filename)
{

}
