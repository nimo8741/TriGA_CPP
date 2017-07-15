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
#include <cstdarg>

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

void nurb::readMsh(std::string filename, int msh_degree)
{
	// set the global mesh degree for later use
	degree = msh_degree;
	
	// calculate the number of nodes which will be present in each triangle
	nodes_in_triangle = ((degree + 2) * (degree + 1) >> 1);
	
	get_bary(degree);


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
	vector<vector<double>>node_vec(nodes,vector<double>(3,1));
	
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
				temp[2] = geo_en - 1;                      // I changed this from the phys_en
				global_edges.push_back(temp);
				i++;
			}
			if (i % 2 == 1) {
				getline(infile, line);   // skip every other line for this section because for some reason, there is always a duplicate line.
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
		current->controlP.resize(nodes_in_triangle);

		// 
		/*

		3
		| \
		|    \
		|       \
		2p+2       2p + 1
		|             \
		|                \
		|        3p+3       \
		|                      \
		|                         \
		|                            \
		|                                \
		3p      3p+1            3p+2         p + 3 
		|                                        \
		|                                            \
		1-------4-----------------------p + 2 ----------2
		

		*/

		vector<vector<double>> current_controlP;
		vector<double> point(3, 1);     // this will set it so that the default weight for the nodes is 1
		point[0] = node_list[node1][0];   // point 1
		point[1] = node_list[node1][1];
		current->controlP[0] = node1;
		current_controlP.push_back(point);


		point[0] = node_list[node2][0];   // point 2
		point[1] = node_list[node2][1];
		current->controlP[1] = node2;
		current_controlP.push_back(point);


		point[0] = node_list[node3][0];   // point 3
		point[1] = node_list[node3][1];
		current->controlP[2] = node3;
		current_controlP.push_back(point);

		vector<double> point1 = node_list[current->controlP[0]]; // have a way to keep track of each of the three vertices of the virtual triangle
		vector<double> point2 = node_list[current->controlP[1]];
		vector<double> point3 = node_list[current->controlP[2]];
		
		
		
		for (int k = 3; k < nodes_in_triangle; k++) {    // start at 3 since I already have the first 3 points
			point[0] = (bary_template[k][0] * point1[0]) + (bary_template[k][1] * point2[0]) + (bary_template[k][2] * point3[0]);
			point[1] = (bary_template[k][0] * point1[1]) + (bary_template[k][1] * point2[1]) + (bary_template[k][2] * point3[1]);
			current_controlP.push_back(point);
		}


		// now go back through and add all of the internal nodes to the node_list since they cannot be repeats
		for (int j = 3 * degree; j < nodes_in_triangle; j++) {
			node_list.push_back(current_controlP[j]);
			current->controlP[j] = int(node_list.size()) - 1;  

		}

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
			else {    // this means that this row does have interior nodes inserted, it will have size degree + 2
				if ((global_edges[j][0] == node1 && global_edges[j][degree] == node2) || (global_edges[j][0] == node2 && global_edges[j][degree] == node1)) {
					// this allows the edge to be saved in either order
					edge1Found = j;

				}
				if ((global_edges[j][0] == node2 && global_edges[j][degree] == node3) || (global_edges[j][0] == node3 && global_edges[j][degree] == node2)) {
					// this allows the edge to be saved in either order
					edge2Found = j;

				}
				if ((global_edges[j][0] == node3 && global_edges[j][degree] == node1) || (global_edges[j][0] == node1 && global_edges[j][degree] == node3)) {
					// this allows the edge to be saved in either order
					edge3Found = j;

				}

			}
			if ((edge1Found != -1) && (edge2Found != -1) && (edge3Found != -1)) {
				break;   // speed things up if they are all found
			}
		}


		///////////////////////////  EDGE 1   /////////////////////////////////
		if (edge1Found == -1) {   // the -1 means that the edge was not found
			// first add the internal nodes to the node_list
			vector <int> temp(degree + 2, 0);
			temp[0] = node1;

			temp[degree] = node2;
			temp[degree + 1] = phys_en;


			for (int k = 3; k < degree + 2; k++) {
				node_list.push_back(current_controlP[k]);
				current->controlP[k] = int(node_list.size()) - 1;
				temp[k - 2] = current->controlP[k];
			}

			// now add the indexes of these nodes into the global edges and global sides fields

			global_edges.push_back(temp); // add this edge to the global list
			current->global_side.push_back(int(global_edges.size()) - 1);   // update the current element's edge parameter with the index of the global edge

		}
		else {
			current->global_side.push_back(edge1Found); // update the global side variable with the index from the global_edge data
			if (global_edges[edge1Found].size() == 3) {  // this means that this edge does not include the p - 1 interior points and so they must be added
				
				for (int k = 3; k < degree + 2; k++) {
					node_list.push_back(current_controlP[k]);
					global_edges[edge1Found].insert(global_edges[edge1Found].begin() + 1, int(node_list.size()) - 1);
					current->controlP[k] = int(node_list.size()) - 1;
				}
			}
			else {   // the global edge does contain the p - 1 interior points
				if (current_controlP[3] == node_list[global_edges[edge1Found][1]]) {   // this means the edge is in the same direction as the global edge is defined
					for (int k = 1; k < degree; k++) {
						current->controlP[k + 2] = global_edges[edge1Found][k];
					}

				}
				else {
					for (int k = 1; k < degree; k++) {
						current->controlP[k + 2] = global_edges[edge1Found][degree - k];
					}
				}

			}
		}


		//////////////////////////   EDGE 2    /////////////////////////////////
		if (edge2Found == -1) {
			// first add the internal nodes to the node_list
			vector <int> temp(degree + 2, 0);
			temp[0] = node2;

			temp[degree] = node3;
			temp[degree + 1] = phys_en;


			for (int k = degree + 2; k < 2*degree + 1; k++) {
				node_list.push_back(current_controlP[k]);
				current->controlP[k] = int(node_list.size()) - 1;
				temp[k - degree - 1] = current->controlP[k];
			}
			global_edges.push_back(temp); // add this edge to the global list
			current->global_side.push_back(int(global_edges.size()) - 1);   // update the current element's edge parameter with the index of the global edge

		}
		else {
			current->global_side.push_back(edge2Found);
			if (global_edges[edge2Found].size() == 3) {  // this means that this edge does not include the 2 interior points and so they must be added

				for (int k = degree + 2; k < 2*degree + 1; k++) {
					node_list.push_back(current_controlP[k]);
					global_edges[edge2Found].insert(global_edges[edge2Found].begin() + 1, int(node_list.size()) - 1);
					current->controlP[k] = int(node_list.size()) - 1;
				}
			}
			else {
				if (current_controlP[degree + 2] == node_list[global_edges[edge2Found][1]]) {   // this means the edge is in the same direction as the global edge is defined
					for (int k = 1; k < degree; k++) {
						current->controlP[degree + 1 + k] = global_edges[edge2Found][k];
					}

				}
				else {
					for (int k = 1; k < degree; k++) {
						current->controlP[degree + 1 + k] = global_edges[edge2Found][degree - k];
					}
				}
			}
		}

		//////////////////////////   EDGE 3    /////////////////////////////////
		if (edge3Found == -1) {
			// first add the internal nodes to the node_list
			vector <int> temp(degree + 2, 0);
			temp[0] = node3;

			temp[degree] = node1;
			temp[degree + 1] = phys_en;


			for (int k = 2*degree + 1; k < 3 * degree; k++) {
				node_list.push_back(current_controlP[k]);
				current->controlP[k] = int(node_list.size()) - 1;
				temp[k - degree - degree] = current->controlP[k];   // subtraction is faster than doing the multiplication of 2* degree
			}
			global_edges.push_back(temp); // add this edge to the global list
			current->global_side.push_back(int(global_edges.size()) - 1);   // update the current element's edge parameter with the index of the global edge

		}
		else {
			current->global_side.push_back(edge3Found);
			if (global_edges[edge3Found].size() == 3) {  // this means that this edge does not include the 2 interior points and so they must be added

				for (int k = 2*degree + 1; k < 3 * degree; k++) {
					node_list.push_back(current_controlP[k]);
					global_edges[edge3Found].insert(global_edges[edge3Found].begin() + 1, int(node_list.size()) - 1);
					current->controlP[k] = int(node_list.size()) - 1;
				}

			}
			else {
				if (current_controlP[2*degree + 1] == node_list[global_edges[edge3Found][1]]) {   // this means the edge is in the same direction as the global edge is defined
					for (int k = 1; k < degree; k++) {
						current->controlP[2*degree + k] = global_edges[edge3Found][k];
					}
				}
				else {
					for (int k = 1; k < degree; k++) {
						current->controlP[2 * degree + k] = global_edges[edge3Found][degree - k];
					}
				}
			}
		}

		// now add current to the list of triangular elements
		triangles.push_back(current);
		getline(infile, line);
	}

}



vector<double> nurb::find_intersect(vector<double> line1p1, vector<double> line1p2, vector<double> line2p1, vector<double> line2p2)
{
	vector<double> intersect(2, 0);
	double slope1, slope2;
	if (line1p1[0] == line1p2[0]) { // there is an infinite slope since the two points lie on top of each other
		intersect[0] = line1p1[0];
		// now solve for the y coor
		// the other one is gaurenteed to have a define slope since the two lines are guarenteed to intersect
		slope2 = (line2p2[1] - line2p1[1]) / (line2p2[0] - line2p1[0]);
		// now using the point-slope form of the equation for a line to find the y coor of the intersection
		intersect[1] = slope2 * (intersect[0] - line2p1[0]) + line2p1[1];
	}
	else if (line2p1[0] == line2p2[0]) {  // line 2 is vertical so I need to rely on line 1
		intersect[0] = line2p1[0];
		slope1 = (line1p2[1] - line1p1[1]) / (line1p2[0] - line1p1[0]);
		intersect[1] = slope1 * (intersect[0] - line1p1[0]) + line1p1[1];
	}
	else {  // now regular case where neither of the slopes are undefined
		slope1 = (-1) * (line1p2[1] - line1p1[1]) / (line1p2[0] - line1p1[0]);  // I have the negative sign because when solving for intersection, I always need the negative of the slope
		slope2 = (-1) * (line2p2[1] - line2p1[1]) / (line2p2[0] - line2p1[0]);
		// now assemble the matrices for eigen to solve
		Matrix2d A;
		Vector2d B;
		A << slope1, 1, slope2, 1;
		B << line1p1[1] + slope1 * line1p1[0], line2p1[1] + slope2 * line2p1[0];
		Vector2d x = A.colPivHouseholderQr().solve(B);    // solve the Ax = B problem
		intersect[0] = x(0);
		intersect[1] = x(1);

	}


	return intersect;
}


void nurb::get_bary(int degree)
{
	bary_template.resize(nodes_in_triangle);
	vector<double> bary_row(3, 0);
	for (int i = 0; i <= degree; i++) {
		for (int j = 0; j <= (degree - i); j++) {
			double s = double(i) / double(degree);
			double r = double(j) / double(degree);
			double t = 1 - s - r;

			// now perform conditioning on s,t, and r
			s  = (s * degree) + 0.5;
			int i_s = (int) s;
			s = (double)i_s;
			s = s / double(degree);

			t = (t * degree) + 0.5;
			int i_t = (int)t;
			t = (double)i_t;
			t = t / double(degree);

			r = (r * degree) + 0.5;
			int i_r = (int)r;
			r = (double)i_r;
			r = r / double(degree);



			bary_row[0] = t;
			bary_row[1] = r;
			bary_row[2] = s;
			int index = ij_to_index(i, j, degree, bary_row);
			bary_template[index] = bary_row;
		}
	}
}

int nurb::ij_to_index(int i, int j, int degree, vector<double> bary)
{
	// handle the case if it is the bounding vertices
	if (bary[0] == 1)
		return 0;
	else if (bary[1] == 1)
		return 1;
	else if (bary[2] == 1)
		return 2;

	// handle the case if it is on the bounding edges

	if (bary[0] == 0) { // this means it is on side 1
		return degree + i + 1;
	}
	else if (bary[1] == 0) {  // this means it is on side 2
		return (degree << 1) + (degree - i);
	}
	else if (bary[2] == 0) {     // this means it is on side 0
		return j + 2;
	}

	// if it made it this far, it is an internal node
	int nodes_checked = 3 * degree - 1;
	while (nodes_checked < nodes_in_triangle) {
		// set up the internal triangle by translating and changing the degree
		i--;
		j--;
		int virtual_degree = degree - 3;
		int nodes_in_level = ((virtual_degree + 1) * (virtual_degree + 2)) >> 1;

		// check for the vertices first
		if (i == 0 && j == 0)
			return nodes_checked + 1;
		else if (i == 0 && j == virtual_degree)
			return nodes_checked + 2;
		else if (i == virtual_degree && j == 0)
			return nodes_checked + 3;

		// now check for sides
		if (i == 0) {   // this means it is on side 0
			return nodes_checked + j + 3;
		}
		else if (j == 0) {   // this means it is on side 2
			return nodes_checked + (virtual_degree - i) + (virtual_degree << 1) + 1;
		}
		else if (i + j == virtual_degree) {   // this means it is on side 1
			return nodes_checked + virtual_degree + i + 2;
		}
		nodes_checked += 3 * virtual_degree;
	}
	return (-1);


}