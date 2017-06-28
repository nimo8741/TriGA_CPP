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
				geo_file << line << endl;

				sprintf_s(line, size, "Physical Line(\"Knot Vector Section %d\") = {%d, %d};", l_count, j + p_start, p_start);
				geo_file << line << endl;

			}
			else {
				sprintf_s(line, size, "Line(%d) = {%d, %d};", l_count, j + p_start, j + 1 + p_start);
				geo_file << line << endl;
				sprintf_s(line, size, "Physical Line(\"Knot Vector Section %d\") = {%d, %d};", l_count, j + p_start, j + 1 + p_start);
				geo_file << line << endl;

			}

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

		/*vector<vector<int>> boundary_group (1, vector<int>(1,0));   // this initializes a dynamic 2D vector with a single 1x1 entry contains a 0 for element 0

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
			string phy_line = "Physical Line(\"Knot Vector Section ";
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
		*/

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
void nurb::call_gmsh(string filename)

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
	// setup the arguments
	
	string args = filename + ".geo -1 -2 -o ";
	args += filename;
	args += ".msh";

	const int size = 200;
	char runline[size];
	sprintf_s(runline, size, "gmsh.exe %s", args.c_str());

	system((char *) runline);   // this line runs gmsh with default settings

}


/**********************************************************************************
Function prototype:
void nurb::readMsh(string filename)

Function description:
This function reads the previously generated .msh file and incorpoorates it into the
program so it can be used for smoothing and degree elevation.

Precondition:
This function uses the filename used in the creation of the .geo file

Postcondition:
This function creates the instances of the triangles within the linear mesh, as well
as a global edge vector which describes which nodes make up which global edges.

This function also insert control points into each mesh triangle so that there is
a total of 10 points.
 
**********************************************************************************/

void nurb::readMsh(std::string filename)
{
	ifstream infile;
	filename += ".msh";
	infile.open(filename);
	string line;       // this will be a reused variable that will hold one line at a time
	double cur_num;     // this will be a reused variable that while help with reading in the numbers

	// Get rid of the initial header lines
	getline(infile, line);    //MeshFormat line
	getline(infile, line);    // format specification line, doesn't really matter
	getline(infile, line);
	getline(infile, line);

	getline(infile, line);
	phy_groups = stoi(line);
	// now loop through all of the physical group data info, I don't really care about this since it is already known
	for (int i = 0; i < phy_groups; i++) {
		getline(infile, line);
	}

	getline(infile, line);    // $EndPhysicalNames
	getline(infile, line);    // $Nodes


	  //*********************************************************//
	  //***********  Loop through all the Nodes   ***************//
	  //*********************************************************//

	getline(infile, line);
	int nodes = stoi(line);
	vector<vector<double>>node_vec(nodes,vector<double>(2,0));
	
	// this will initialize a 2D vector that contains enough rows to hold all of the nodes and 2 columns for the X and Y position
	// this data will be used later when constructing stucts which contain the data for each of the elements. The format of the .msh
	// file is what forces be to define this temporary vector array.


	for (int i = 0; i < nodes; i++) {
		getline(infile, line);    // grab the next line
		int k = 0;
		int l;
		int j = 0;
		while (j < 3) {
			// loop through all of the white space
			while (line[k] == ' ') {
				k = k++;
			}

			// Loop through all the numbers
			l = 0;
			while (line[k + l] != ' ') {
				l = l++;
			}
			string temp_num = line.substr(k, l);
			cur_num = atof(temp_num.c_str());
			switch (j) {
			case 0:
				// do nothing since this is just the node number
				break;
			case 1:
				node_vec[i][0] = cur_num;
				break;
			case 2:
				node_vec[i][1] = cur_num;
				break;
			}
			j++;
			k = k + l;
		}

	}
	node_list = node_vec;

	//*********************************************************//
	//***********  Loop through all the Elements    ***********//
	//*********************************************************//

	getline(infile, line);      // $EndNodes
	getline(infile, line);      // $Elements
	getline(infile, line);
	int num_elem = stoi(line);

	/* The format for the element lines are as follows

	<Element Number>, <Element Type Code>, <number of "Tags">, first the physical entity the element belongs in, second the elementary geometrical entity, <node number list>
	Element type codes
	1: 2 node line
	2: 3 node triangle
	*/
	int type = 1;
	for (int i = 0; i < num_elem; i++) { 
		if (i > 71)
			bool fuck = true;
		// first I need to skip through all of the elements that are 2 node lines as they are not important
		int index, num_tag, phys_en, geo_en, node1, node2, node3;
		while (type == 1) {
			getline(infile, line);
			stringstream ss(line);
			ss >> index >> type >> num_tag >> phys_en >> geo_en >> node1 >> node2;
			// now since all of these line element are distinct, I can just add all of them to the global_edge vector
			// since if they shared an edge, they would have to be the same element.
			if (type == 1) {
				vector<int> temp(3, 0);
				temp[0] = node1 - 1;
				temp[1] = node2 - 1;
				temp[2] = phys_en;
				global_edges.push_back(temp);
				i++;
			}
		}
		// now I am left with only the triangles
		stringstream ss(line);
		ss >> index >> type >> num_tag >> phys_en >> geo_en >> node1 >> node2 >> node3;

		node1--;
		node2--;   // I am subracting 1 for all of these so that they match the (starting at 0) indexes of node_list
		node3--;

		// now put this data into the tri_elem structs
		Tri_elem *current = new Tri_elem;
		vector<int> IEN_row(10, 0);   // this will comprise the row in the IEN array, so I will make it while I go through this next section of code and then add it to the full IEN

		// first deal with the 10 control points, also I will do some of the loop unrolling myself here for efficiency the control points are laid out like
		/*

		3
		| \
		|  \
		8   7
		|    \
		|     \
		9  10  6
		|       \
		|        \
		1--4---5--2
		

		*/


		vector<double> point(2, 0);
		point[0] = node_list[node1][0];   // point 1
		point[1] = node_list[node1][1];
		current->controlP.push_back(point);
		IEN_row[0] = node1;

		point[0] = node_list[node2][0];   // point 2
		point[1] = node_list[node2][1];
		current->controlP.push_back(point);
		IEN_row[1] = node2;

		point[0] = node_list[node3][0];   // point 3
		point[1] = node_list[node3][1];
		current->controlP.push_back(point);
		IEN_row[2] = node3;

		// now figure out the interpolation points there are found at each third down the length of the edge
		point[0] = (2.0 / 3.0) * current->controlP[0][0] + (1.0 / 3.0) * current->controlP[1][0];
		point[1] = (2.0 / 3.0) * current->controlP[0][1] + (1.0 / 3.0) * current->controlP[1][1];    // point 4
		current->controlP.push_back(point);

		point[0] = (1.0 / 3.0) * current->controlP[0][0] + (2.0 / 3.0) * current->controlP[1][0];
		point[1] = (1.0 / 3.0) * current->controlP[0][1] + (2.0 / 3.0) * current->controlP[1][1];    // point 5
		current->controlP.push_back(point);
		/////////////////////////////////////////////////////////////////////////////////////////////

		point[0] = (2.0 / 3.0) * current->controlP[1][0] + (1.0 / 3.0) * current->controlP[2][0];
		point[1] = (2.0 / 3.0) * current->controlP[1][1] + (1.0 / 3.0) * current->controlP[2][1];    // point 6
		current->controlP.push_back(point);

		point[0] = (1.0 / 3.0) * current->controlP[1][0] + (2.0 / 3.0) * current->controlP[2][0];
		point[1] = (1.0 / 3.0) * current->controlP[1][1] + (2.0 / 3.0) * current->controlP[2][1];    // point 7
		current->controlP.push_back(point);
		/////////////////////////////////////////////////////////////////////////////////////////////

		point[0] = (2.0 / 3.0) * current->controlP[2][0] + (1.0 / 3.0) * current->controlP[0][0];
		point[1] = (2.0 / 3.0) * current->controlP[2][1] + (1.0 / 3.0) * current->controlP[0][1];    // point 8
		current->controlP.push_back(point);

		point[0] = (1.0 / 3.0) * current->controlP[2][0] + (2.0 / 3.0) * current->controlP[0][0];
		point[1] = (1.0 / 3.0) * current->controlP[2][1] + (2.0 / 3.0) * current->controlP[0][1];    // point 9
		current->controlP.push_back(point);
		/////////////////////////////////////////////////////////////////////////////////////////////

		point[0] = (1.0 / 2.0) * (current->controlP[5][0] + current->controlP[8][0]);
		point[1] = (1.0 / 2.0) * (current->controlP[5][1] + current->controlP[8][1]);  // point 10, this is the midpoint between point 6 and 9
		current->controlP.push_back(point);
		// since this point 10 is guarenteed to be unique, I can add it to the node_list
		node_list.push_back(point);
		IEN_row[9] = node_list.size() - 1;

		//****************************************************************//
		//************* now figure out the global edge indexes************//
		//****************************************************************//


		// first figure out if the edge already exists and if it doesn't add it to the global list
		int edge1Found = -1;
		int edge2Found = -1;
		int edge3Found = -1;

		for (unsigned int j = 0; j < global_edges.size(); j++) {
			if (global_edges[j].size() == 3) {   // this means that this row does not have the interior nodes inserted
				if ((global_edges[j][0] == node1 && global_edges[j][1] == node2) || (global_edges[j][0] == node2 && global_edges[j][1] == node1)) {
					// this allows the edge to be saved in either order
					edge1Found = j;

				}
				if ((global_edges[j][0] == node2 && global_edges[j][1] == node3) || (global_edges[j][0] == node3 && global_edges[j][1] == node2)) {
					// this allows the edge to be saved in either order
					edge2Found = j;

				}
				if ((global_edges[j][0] == node3 && global_edges[j][1] == node1) || (global_edges[j][0] == node1 && global_edges[j][1] == node3)) {
					// this allows the edge to be saved in either order
					edge3Found = j;

				}
			}
			else {    // this means that this row does have interior nodes inserted, it will have size 5
				if ((global_edges[j][0] == node1 && global_edges[j][3] == node2) || (global_edges[j][0] == node2 && global_edges[j][3] == node1)) {
					// this allows the edge to be saved in either order
					edge1Found = j;

				}
				if ((global_edges[j][0] == node2 && global_edges[j][3] == node3) || (global_edges[j][0] == node3 && global_edges[j][3] == node2)) {
					// this allows the edge to be saved in either order
					edge2Found = j;

				}
				if ((global_edges[j][0] == node3 && global_edges[j][3] == node1) || (global_edges[j][0] == node1 && global_edges[j][3] == node3)) {
					// this allows the edge to be saved in either order
					edge3Found = j;

				}

			}
			if ((edge1Found != -1) && (edge2Found != -1) && (edge3Found != -1)) {
				break;   // speed things up if they are all found
			}
		}


		///////////////////////////  EDGE 1   /////////////////////////////////
		if (edge1Found == -1) {   // the -1 means that the node was not found
			// first add the internal nodes to the node_list
			node_list.push_back(current->controlP[3]);
			IEN_row[3] = node_list.size() - 1;   // point 4

			node_list.push_back(current->controlP[4]);
			IEN_row[4] = node_list.size() - 1;    // point 5

			// now add the indexes of these nodes into the global edges and global sides fields
			vector <int> temp(5, 0);
			temp[0] = node1;
			temp[1] = node_list.size() - 2;
			temp[2] = temp[1] + 1;
			temp[3] = node2;
			temp[4] = phys_en;
			global_edges.push_back(temp); // add this edge to the global list
			current->global_side.push_back(global_edges.size() - 1);   // update the current element's edge parameter with the index of the global edge

		}
		else {
			current->global_side.push_back(edge1Found); // update the global side variable with the index from the global_edge data
			if (global_edges[edge1Found].size() == 3) {  // this means that this edge does not include the 2 interior points and so they must be added
				
				// node 4
				node_list.push_back(current->controlP[3]);
				global_edges[edge1Found].insert(global_edges[edge1Found].begin() + 1, node_list.size() - 1);
				IEN_row[3] = node_list.size() - 1;

				// node 5
				node_list.push_back(current->controlP[4]);
				global_edges[edge1Found].insert(global_edges[edge1Found].begin() + 2, node_list.size() - 1);
				IEN_row[4] = node_list.size() - 1;

			}
			else {
				if (current->controlP[3] == node_list[global_edges[edge1Found][1]]) {   // this means the edge is in the same direction as the global edge is defined
					IEN_row[3] = global_edges[edge1Found][1];   // update with node 4
					IEN_row[4] = global_edges[edge1Found][2];     // I'm not sure if it ever goes into here

				}
				else {
					IEN_row[3] = global_edges[edge1Found][2];  // this means that the edge is the reverse order so I need to switch the assignment
					IEN_row[4] = global_edges[edge1Found][1];
				}

			}
		}




		//////////////////////////   EDGE 2    /////////////////////////////////
		if (edge2Found == -1) {
			// first add the internal nodes to the node_list
			node_list.push_back(current->controlP[5]);
			IEN_row[5] = node_list.size() - 1;   // point 6

			node_list.push_back(current->controlP[6]);
			IEN_row[6] = node_list.size() - 1;   // point 7


			vector <int> temp(5, 0);
			temp[0] = node2;
			temp[1] = node_list.size() - 2;
			temp[2] = temp[1] + 1;
			temp[3] = node3;
			temp[4] = phys_en;
			global_edges.push_back(temp); // add this edge to the global list
			current->global_side.push_back(global_edges.size() - 1);   // update the current element's edge parameter with the index of the global edge

		}
		else {
			current->global_side.push_back(edge2Found);
			if (global_edges[edge2Found].size() == 3) {  // this means that this edge does not include the 2 interior points and so they must be added

				// node 6
				node_list.push_back(current->controlP[5]);
				global_edges[edge2Found].insert(global_edges[edge2Found].begin() + 1, node_list.size() - 1);
				IEN_row[5] = node_list.size() - 1;

				// node 7
				node_list.push_back(current->controlP[6]);
				global_edges[edge2Found].insert(global_edges[edge2Found].begin() + 2, node_list.size() - 1);
				IEN_row[6] = node_list.size() - 1;
			}
			else {
				if (current->controlP[5] == node_list[global_edges[edge2Found][1]]) {   // this means the edge is in the same direction as the global edge is defined
					IEN_row[5] = global_edges[edge2Found][1];   // update with node 6
					IEN_row[6] = global_edges[edge2Found][2];    // I'm not sure if it ever goes into here

				}
				else {
					IEN_row[5] = global_edges[edge2Found][2];  // this means that the edge is the reverse order so I need to switch the assignment
					IEN_row[6] = global_edges[edge2Found][1];
				}
			}
		}

		//////////////////////////   EDGE 3    /////////////////////////////////
		if (edge3Found == -1) {
			// first add the internal nodes to the node_list
			node_list.push_back(current->controlP[7]);
			IEN_row[7] = node_list.size() - 1;   // point 8

			node_list.push_back(current->controlP[8]);
			IEN_row[8] = node_list.size() - 1;   // point 9


			vector <int> temp(5, 0);
			temp[0] = node3;
			temp[1] = node_list.size() - 2;
			temp[2] = temp[1] + 1;
			temp[3] = node1;
			temp[4] = phys_en;
			global_edges.push_back(temp); // add this edge to the global list
			current->global_side.push_back(global_edges.size() - 1);   // update the current element's edge parameter with the index of the global edge

		}
		else {
			current->global_side.push_back(edge3Found);
			if (global_edges[edge3Found].size() == 3) {  // this means that this edge does not include the 2 interior points and so they must be added

				 // node 8
				node_list.push_back(current->controlP[7]);
				global_edges[edge3Found].insert(global_edges[edge3Found].begin() + 1, node_list.size() - 1);

				// node 9
				node_list.push_back(current->controlP[8]);
				global_edges[edge3Found].insert(global_edges[edge3Found].begin() + 2, node_list.size() - 1);

			}
			else {
				if (current->controlP[7] == node_list[global_edges[edge3Found][1]]) {   // this means the edge is in the same direction as the global edge is defined
					IEN_row[7] = global_edges[edge3Found][1];   // update with node 8
					IEN_row[8] = global_edges[edge3Found][2];      // I'm not sure if it ever goes into here

				}
				else {
					IEN_row[7] = global_edges[edge3Found][2];  // this means that the edge is the reverse order so I need to switch the assignment
					IEN_row[8] = global_edges[edge3Found][1];
				}
			}
		}

		// now add current to the list of triangular elements
		triangles.push_back(current);
		Tri_IEN.push_back(IEN_row);
		getline(infile, line);
	}

}

