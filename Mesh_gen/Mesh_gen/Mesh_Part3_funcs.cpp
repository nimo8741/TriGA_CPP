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
#include <random>
#include <numeric>
#include <list>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void nurb::smoothMesh(int mesh_degree)
{
	degree = mesh_degree;
	create_side_nodes();
	organize_boundary();
	split_and_extract();
	assign_boundary_points();
	bool fuck = true;

	//boundary_weights(degree);
	// now I need to solve for the weights of all of the points]
	//smooth_weights(degree);

	// now I need to solve the linear elasticity
	//LE2D(degree);
}

void nurb::create_side_nodes()
{
	vector<int>node_side1(degree + 1);
	node_side1[0] = 0;
	node_side1[degree] = 1;
	for (int j = 1; j < degree; j++) {
		node_side1[j] = j + 2;
	}

	vector<int>node_side2(degree + 1);
	node_side2[0] = 1;
	node_side2[degree] = 2;
	for (int j = 1; j < degree; j++) {
		node_side2[j] = j + degree + 1;
	}

	vector<int>node_side3(degree + 1);
	node_side3[0] = 2;
	node_side3[degree] = 0;
	for (int j = 1; j < degree; j++) {
		node_side3[j] = j + (degree << 1);
	}
	node_side_index.push_back(node_side1);
	node_side_index.push_back(node_side2);
	node_side_index.push_back(node_side3);
}

void nurb::organize_boundary()
{
	for (unsigned int i = 0; i < triangles.size(); i++) {  //  loop through all of the triangles
		for (int side = 0; side < 3; side++) {
			if (global_edges[triangles[i]->global_side[side]][degree + 1] != phy_groups) {  // if it passes this, it is a boundary edge
				// now I need to loop through all of the nodes that would be on the boundary, there are degree + 1 of them

				/////////////////////////////////////////////////////////////////////////////////
				//Create a list of all of the node indexes where the node's are on the boundary
				/////////////////////////////////////////////////////////////////////////////////

				for (int edge_node = 0; edge_node < degree + 1; edge_node++) {
					bool found = false;
					for (unsigned int j = 0; j < bNodes.size(); j++) {
						if (bNodes[j] == triangles[i]->controlP[node_side_index[side][edge_node]]) {
							found = true;
							break;
						}
					}
					// if the node was not found within bNodes, add it to the list
					if (!found)
						bNodes.push_back(triangles[i]->controlP[node_side_index[side][edge_node]]);

				}

				////////////////////////////////////////////////////////////////////////////////////////
				//now determine which of the element section the current boundary triangle is part of
				/////////////////////////////////////////////////////////////////////////////////////////

				unsigned int KV_section = global_edges[triangles[i]->global_side[side]][degree + 1];   // KV_section starts at 1
				unsigned int cur_nurb = 0;   // this will be the identifier for the current nurbs curve.
				bool found = false;
				while (found == false) {
					if (cur_nurb == 1)
						bool fuck = true;
					// in the way the mesh is defined, there is a straight line connecting all of the xi_real_loc point for the NURBS curve

					unsigned xi_size = Elem_list[cur_nurb]->xi_evals.size();
					if (KV_section < Elem_list[cur_nurb]->n_el * (xi_size - 1)) { // this means the KV_section lies in the current NURBS curve
						unsigned int element = KV_section / (xi_size - 1); // this is the identifier within the bezier element for that NURBS curve. Starts at 0
						unsigned int elem_section = KV_section % (xi_size - 1);    // this is the identifier for the section within that bezier element.  Starts at 0

						vector<unsigned int> tri_NURB_row(5, 0);
						tri_NURB_row[0] = i;
						tri_NURB_row[1] = cur_nurb;
						tri_NURB_row[2] = element;
						tri_NURB_row[3] = elem_section;
						tri_NURB_row[4] = side;
						tri_NURB_elem_section_side.push_back(tri_NURB_row);

						found = true;
					}
					else {
						KV_section -= Elem_list[cur_nurb]->n_el * (xi_size - 1);
						cur_nurb++;
					}
				}
			}
		}
	}
	unsigned int size = tri_NURB_elem_section_side.size();
	vector<int> work_around(size);
	iota(begin(work_around), end(work_around), 0);
	sort (work_around.begin(), work_around.end(), *this);
	vector<vector<unsigned int>> temp(size);
	for (unsigned int i = 0; i < size; i++) {
		temp[i] = tri_NURB_elem_section_side[work_around[i]];
	}
	tri_NURB_elem_section_side = temp;
}

bool nurb::operator()(int i, int j)
{
	vector<unsigned int> input2 = tri_NURB_elem_section_side[j];
	vector<unsigned int> input1 = tri_NURB_elem_section_side[i];
	if (input2[1] > input1[1])
		return true;
	else if (input2[1] < input1[1])
		return false;
	else {  // this means that the nurbs curves are equal
		if (input2[2] > input1[2])
			return true;
		else if (input2[2] < input1[2])
			return false;
		else{   // this means that the element indexes are the same
			if (input2[3] > input1[3])
				return true;
			else if (input2[3] < input1[3])
				return false;
			else {     // this means that the sections are the same so now I need to order within 

				vector<double> anchor = Elem_list[input1[1]]->elem_Geom[input1[2]].controlP[0];
				vector<double> first_point1 = node_list[triangles[input1[0]]->controlP[node_side_index[input1[4]][0]]];
				vector<double> second_point1 = node_list[triangles[input1[0]]->controlP[node_side_index[input1[4]][degree]]];;

				vector<double> first_point2 = node_list[triangles[input2[0]]->controlP[node_side_index[input2[4]][0]]];;
				vector<double> second_point2 = node_list[triangles[input2[0]]->controlP[node_side_index[input2[4]][degree]]];;


				double dist = sqrt(pow(first_point1[0] - anchor[0], 2) + pow(first_point1[1] - anchor[1], 2));
				double temp = sqrt(pow(second_point1[0] - anchor[0], 2) + pow(second_point1[1] - anchor[1], 2));
				dist = min(dist, temp);

				double dist2 = sqrt(pow(first_point2[0] - anchor[0], 2) + pow(first_point2[1] - anchor[1], 2));
				temp = sqrt(pow(second_point2[0] - anchor[0], 2) + pow(second_point2[1] - anchor[1], 2));
				dist2 = min(dist2, temp);

				if (dist2 > dist)
					return true;
				else      // dist2 < dist
					return false;


			}
		}
	}
}

void nurb::split_and_extract()
{
	// I will need to do this for each of the elements within each of the boundary curves
	for (int cur_nurb = 0; cur_nurb < num_curves; cur_nurb++) {
		for (int cur_elem = 0; cur_elem < Elem_list[cur_nurb]->n_el; cur_elem++) {
			vector<double> xi_to_add = determine_elem_split(cur_nurb, cur_elem);
			int num_to_increment = xi_to_add.size();
			// now I need to add enough knots at these locations so that there is c(0) continuity
			int nurb_degree = Elem_list[cur_nurb]->p;
			for (unsigned int i = 0; i < xi_to_add.size(); i++) {
				vector<double> repeats(nurb_degree - 1, xi_to_add[i]);
				xi_to_add.insert(xi_to_add.begin() + i, repeats.begin(), repeats.end());
				i += repeats.size();
			}
			// now call curve_refine to refine, split, and elevate the element

			curve_refine(Elem_list[cur_nurb], cur_elem, xi_to_add, false);

			// now need to update tri_NURB_elem_section
			// first I need to find the first element that needs to be changed
			bool found = false;
			unsigned int i = 0;
			while (!found) {
				if (tri_NURB_elem_section_side[i][1] == cur_nurb && tri_NURB_elem_section_side[i][2] == cur_elem)
					found = true;
				else
					i++;
			}
			// now I have the index of the first section within the list since it was previously orgainized
			// I need to update all of the sections for the rest of the element, making them their own elements and setting the section to 0
			i++;
			int elem_index = cur_elem;
			while (i < tri_NURB_elem_section_side.size() && tri_NURB_elem_section_side[i][1] == cur_nurb && tri_NURB_elem_section_side[i][2] == cur_elem) {
				elem_index++;
				tri_NURB_elem_section_side[i][2] = elem_index;
				tri_NURB_elem_section_side[i][3] = 0;
				i++;
			}
			// now I need to loop through all of the remaining elements in the nurb, increasing their element number by elem_index
			while (i < tri_NURB_elem_section_side.size() && tri_NURB_elem_section_side[i][1] == cur_nurb) {
				tri_NURB_elem_section_side[i][2] += num_to_increment;
				i++;
			}


			cur_elem += num_to_increment;
		}
	}

}


vector<double> nurb::determine_elem_split(int cur_nurb, int cur_elem)
{
	// tri_NURB_elem_section_side is already organized so that first everything is of the same curve, then the same element, then increasing section position

	int i = 0;
	while (tri_NURB_elem_section_side[i][1] != cur_nurb)
		i++;
	while (tri_NURB_elem_section_side[i][2] != cur_elem)
		i++;
	int i_save = i;
	double total_dist = 0.0;
	vector<double> split_list;
	unsigned int max_i = tri_NURB_elem_section_side.size();
	while (unsigned (i) < max_i && tri_NURB_elem_section_side[i][2] == cur_elem && tri_NURB_elem_section_side[i][1] == cur_nurb) {
		vector<double> first_point = node_list[triangles[tri_NURB_elem_section_side[i][0]]->controlP[node_side_index[tri_NURB_elem_section_side[i][4]][0]]];
		vector<double> second_point = node_list[triangles[tri_NURB_elem_section_side[i][0]]->controlP[node_side_index[tri_NURB_elem_section_side[i][4]][degree]]];;
		double dist = sqrt(pow(first_point[0] - second_point[0], 2) + pow(first_point[1] - second_point[1], 2));
		split_list.push_back(dist);
		total_dist += dist;
		i++;
	}
	for (unsigned int j = 0; j < split_list.size() - 1; j++) {
		if (j == 0) {
			split_list[j] = split_list[j] / total_dist;
		}
		else {
			split_list[j] = (split_list[j] / total_dist) + split_list[j - 1];
		}
	}
	split_list.pop_back();
	return split_list;
}

void nurb::curve_refine(Bezier_handle * Bez, int cur_elem, vector<double> xi_to_add, bool part1)
{
	// Step 1: create a knot vector for this element.  Since it is Bezier, it is easy to compute
	int p = Bez->p;
	unsigned int n = Bez->elem_Geom[cur_elem].controlP.size();
	vector<double> KV;
	for (int i = 0; i < (p + 1) << 1; i++) {
		KV.push_back(double(i / (p + 1)));
	}

	// Step 2: Project the control points
	for (unsigned int i = 0; i < n; i++) {
		Bez->elem_Geom[cur_elem].controlP[i][0] *= Bez->elem_Geom[cur_elem].weight[i];
		Bez->elem_Geom[cur_elem].controlP[i][1] *= Bez->elem_Geom[cur_elem].weight[i];
		if (Bez->elem_Geom[cur_elem].controlP[i].size() != 3)  // this will either be 2 or 3
			Bez->elem_Geom[cur_elem].controlP[i].push_back(Bez->elem_Geom[cur_elem].weight[i]);
	}

	vector<vector<double>> P = Bez->elem_Geom[cur_elem].controlP;
	vector<vector<double>> Q = P;
	vector<double> KV_temp = KV;
	// Step 3: Knot Insertion algorithm
	for (unsigned int i = 0; i < xi_to_add.size(); i++) {
		double add_cur = xi_to_add[i];
		KV_temp.push_back(add_cur);
		sort(KV_temp.begin(), KV_temp.end());
		ptrdiff_t index = find(KV_temp.begin(), KV_temp.end(), add_cur) - KV_temp.begin();
		unsigned int k = index - 1;
		unsigned int m = n + 1;
		double a;
		for (unsigned int j = 0; j < m; j++) {
			// figure out a
			if (j < k - p + 1)
				a = 1;
			else if (j >= k - p && j < k + 1)
				a = (add_cur - KV[j]) / (KV[j + p] - KV[j]);
			else
				a = 0;

			// Figure out Q_j
			if (!j) {
				Q[j][0] = P[j][0];
				Q[j][1] = P[j][1];
				Q[j][2] = P[j][2];
			}
			else if (j && j < n) {
				Q[j][0] = a*P[j][0] + (1 - a)*P[j - 1][0];
				Q[j][1] = a*P[j][1] + (1 - a)*P[j - 1][1];
				Q[j][2] = a*P[j][2] + (1 - a)*P[j - 1][2];


			}
			else {  // this is where Q gains an additional knot
				vector<double> temp(3, 0);
				temp[0] = P[n - 1][0];
				temp[1] = P[n - 1][1];
				temp[2] = P[n - 1][2];
				Q.push_back(temp);
				n++;  // increment the number of control points

			}
		}
		P = Q;
		KV = KV_temp;
	}
	Bez->elem_Geom[cur_elem].controlP = Q;

	// now I need to seperate the 1 element into multiple
	
	unsigned int i;
	for (i = 0; i <= (xi_to_add.size() / Bez->p); i++) {

		// first i need to erase p elements from the front of the current one
		if (i){
			for (int j = 0; j < (Bez->p); j++) {
				Bez->elem_Geom[cur_elem + i].controlP.erase(Bez->elem_Geom[cur_elem + i].controlP.begin());
			}
		}

		if (Bez->elem_Geom[cur_elem + i].controlP.size() == unsigned(Bez->p + 1))
			break;
		// now, I need to make a copy of the element data and put it following the current element
		Bez->elem_Geom.insert(Bez->elem_Geom.begin() + cur_elem + i + 1, Bez->elem_Geom[cur_elem + i]);

		// now, remove all of the nodes past the p+1 that it need from the current one
		while (Bez->elem_Geom[cur_elem + i].controlP.size() > unsigned(Bez->p + 1)) {
			Bez->elem_Geom[cur_elem + i].controlP.pop_back();
		}

		Bez->n_el++;
	}


	if (!part1) {    // degree elevation in here
		for (unsigned int j = 0; j <= i; j++) {   // loop through all of the elements which need elevating
			int k;
			if (j == i && Bez->n_el - 1 - i == cur_elem) {
				while (degree != Bez->p) {   // if this is the last element in the nurbs curve, it will be increasing Bez->p so I need to have a different conditional than below
					elevate_degree(Bez, cur_elem + j);
				}
			}
			else {
				for (k = 0; k < (degree - Bez->p); k++) {
					elevate_degree(Bez, cur_elem + j);
				}
			}
			for (unsigned int k = 0; k < Bez->elem_Geom[cur_elem + j].controlP.size(); k++) {  // unproject the control points

				Bez->elem_Geom[cur_elem + j].controlP[k][0] = Bez->elem_Geom[cur_elem + j].controlP[k][0] / Bez->elem_Geom[cur_elem + j].controlP[k][2];
				Bez->elem_Geom[cur_elem + j].controlP[k][1] = Bez->elem_Geom[cur_elem + j].controlP[k][1] / Bez->elem_Geom[cur_elem + j].controlP[k][2];

				if (k >= Bez->elem_Geom[cur_elem + j].weight.size())
					Bez->elem_Geom[cur_elem + j].weight.push_back(Bez->elem_Geom[cur_elem + j].controlP[k][2]);
				else
					(Bez->elem_Geom[cur_elem + j].weight[k] = Bez->elem_Geom[cur_elem + j].controlP[k][2]);
			}
		}

	}

	if (part1) {
		//  now I need to unproject the control points
		for (unsigned int j = 0; j <= i; j++) {
			for (int k = 0; k <= Bez->p; k++) {

				Bez->elem_Geom[cur_elem + j].controlP[k][0] = Bez->elem_Geom[cur_elem + j].controlP[k][0] / Bez->elem_Geom[cur_elem + j].controlP[k][2];
				Bez->elem_Geom[cur_elem + j].controlP[k][1] = Bez->elem_Geom[cur_elem + j].controlP[k][1] / Bez->elem_Geom[cur_elem + j].controlP[k][2];
				Bez->elem_Geom[cur_elem + j].weight[k] = Bez->elem_Geom[cur_elem + j].controlP[k][2];
			}

			if (j) {
				//// so simply make a copy of the previous Operator
				Bez->Operator.insert(Bez->Operator.begin() + cur_elem + j, Bez->Operator[cur_elem]);

				// Now Update the master IEN array so it reflects the added control points.
				vector<int> insert(Bez->p + 1, 0);  // this is the entry corresponding to the inserted point
				vector<vector<int>> cur_IEN = Master_IEN.back();
				insert[0] = cur_IEN[cur_elem][Bez->p];
				for (int i = 1; i < Bez->p + 1; i++) {
					insert[i] = insert[i - 1] + 1;
				}
				cur_IEN.insert(cur_IEN.begin() + cur_elem + j, insert);

				for (int i = cur_elem + 2; i < Bez->n_el; i++) {
					for (int j = 0; j < Bez->p + 1; j++) {
						cur_IEN[i][j] += Bez->p;
					}

				}
				Master_IEN.back() = cur_IEN;
			}
				// now I need to update the curve and edge lengths of these new elements
				Geom_data *garbage = new Geom_data;    // just need one of these so it will be happy
				vector<double> trash(1, 0);
				evalBez(garbage, Bez, cur_elem + j, trash);
				curve_len(Bez, cur_elem + j);
				delete garbage;
		}
	}

}

void nurb::elevate_degree(Bezier_handle * Bez, int element)
{
	vector<vector<double>> P = Bez->elem_Geom[element].controlP;
	vector<vector<double>> Pe = P;
	int p = Bez->p;
	unsigned int n = P.size();
	Pe.push_back(Pe[n - 1]);

	for (unsigned int i = 1; i < n; i++) {
		Pe[i][0] = (double(i) / double(p + 1))*P[i - 1][0] + (1 - (double(i) / double(p + 1)))*P[i][0];
		Pe[i][1] = (double(i) / double(p + 1))*P[i - 1][1] + (1 - (double(i) / double(p + 1)))*P[i][1];
		Pe[i][2] = (double(i) / double(p + 1))*P[i - 1][2] + (1 - (double(i) / double(p + 1)))*P[i][2];

	}

	Bez->elem_Geom[element].controlP = Pe;
	
	// now figure out what to make the degree
	unsigned int min_p = Bez->elem_Geom[0].controlP.size() - 1;
	for (int i = 1; i < Bez->n_el; i++) {
		if (Bez->elem_Geom[i].controlP.size() - 1 < min_p)
			min_p = Bez->elem_Geom[i].controlP.size() - 1;
	}
	Bez->p = min_p;

}

void nurb::assign_boundary_points()
{
	// loop through tri_NURB_elem_section_side
	for (unsigned int i = 0; i < tri_NURB_elem_section_side.size(); i++) {

		// loop through all of the points
		for (int j = 0; j <= degree; j++){
			// determine the correct value of the control point

			int tri = tri_NURB_elem_section_side[i][0];
			int nurb = tri_NURB_elem_section_side[i][1];
			int elem = tri_NURB_elem_section_side[i][2];
			int side = tri_NURB_elem_section_side[i][4];
			int node = node_side_index[side][j];

			vector<double> actual_p(3, 0);
			actual_p[0] = Elem_list[nurb]->elem_Geom[elem].controlP[j][0];
			actual_p[1] = Elem_list[nurb]->elem_Geom[elem].controlP[j][1];
			actual_p[2] = Elem_list[nurb]->elem_Geom[elem].controlP[j][2];
			
			// determine which node index this is
			int node_list_index = triangles[tri]->controlP[node];

			// update the list
			node_list[node_list_index] = actual_p;
		}
	}
}

/*
void nurb::boundary_weights(int degree)
{
	// Instead of thinking in triangles, I can instead again visualize the boundary as a NURBS curve
	// First, I need to get all of the triangular element which are
	//		1) on the same boundary
	//		2) exist in the same bezier element
	//		3) sort these in order from element section 0 to the end
	//
	// Second, I will need to figure out the parametric xi values for each of the nodes on this curve
	//
	// Third, project the coordinates of the nodes I know the weights for
	//
	// Fourth, perfrom curve refinement and most likely degree elevation to reach the final curve
	//
	// Fifth "un-project" the coordinates to get the weights of all of the intermittent points

	// create a lookup table which will be used later in the function
	// create vector which will tell the follow loop which nodes to grab, this depends on the side on degree of the triangle
	vector<int> node_side1(degree + 1, 0);
	node_side1[0] = 0;
	node_side1[degree] = 1;
	for (int j = 1; j < degree; j++) {
		node_side1[j] = j + 2;
	}

	vector<int> node_side2(degree + 1, 0);
	node_side2[0] = 1;
	node_side2[degree] = 2;
	for (int j = 1; j < degree; j++) {
		node_side2[j] = j + degree + 1;
	}

	vector<int> node_side3(degree + 1, 0);
	node_side3[0] = 2;
	node_side3[degree] = 0;
	for (int j = 1; j < degree; j++) {
		node_side3[j] = j + (degree << 1);
	}


	/////////////////////   STEP 1   ////////////////////////////
	for (int cur_nurb = 0; cur_nurb < num_curves; cur_nurb++) {
		for (int cur_elem = 0; cur_elem < Elem_list[cur_nurb]->n_el; cur_elem++) {
			vector<unsigned int> same_elem_handle;
			for (unsigned int i = 0; i < tri_NURB_elem_section_side.size(); i++) {
				if ((tri_NURB_elem_section_side[i][1] == cur_nurb) && (tri_NURB_elem_section_side[i][2] == cur_elem)) {
					if (same_elem_handle.size() == 0)   // the size is 0 so I can add it in without care
						same_elem_handle.push_back(i);
					else if (same_elem_handle.back() < tri_NURB_elem_section_side[i][3])  // the list is in good order so I can slap the new entry onto the back
						same_elem_handle.push_back(i);
					else {     // the list will not be in order so I need to added the new element such that it will be in order
						unsigned int j;
						for (j = 0; j < same_elem_handle.size(); j++) {
							// if I find an element section in the current element list which is larger than the element section I am examining (at index i) then insert accordingly
							// so same_elem_handle is in the appropriate order by element section, smallest to largest
							if (tri_NURB_elem_section_side[same_elem_handle[j]][3] > tri_NURB_elem_section_side[i][3])
								break;
						}
						// now insert at index j
						same_elem_handle.insert(same_elem_handle.begin() + j, i);
					}
				}
			}
			/////////////////////   STEP 2   ////////////////////////////
			// I know that there are degree+1 point along each edge of the triangle and there is a shared node between each adjacent triangle's edges
			// Also I know that all of these points are equispaced down each of the respective sections
			// I'm not sure if each of the sections within the element are of equivalent size so I will assume that they are not

			vector<double> xi_list(1, 0.0);
			// first point of the first node is automatically a zero so I don't need to worry about it
			vector<double> all_section_lens;
			double element_len = 0.0;

			// figure out the length of the section, I will do this by linearally connecting the degree + 1 points
			for (unsigned int section = 0; section < same_elem_handle.size(); section++) {
				double section_len = 0.0;
				vector<double> point1(2, 0.0);
				vector<double> point2(2, 0.0);

				int tri_index = tri_NURB_elem_section_side[same_elem_handle[section]][0];
				vector<int> nodes_to_use;
				if (tri_NURB_elem_section_side[same_elem_handle[section]][4] == 0) {   // this is side 1
					nodes_to_use = node_side1;
				}
				else if (tri_NURB_elem_section_side[same_elem_handle[section]][4] == 1) {    // this is side 2
					nodes_to_use = node_side2;
				}
				else {     // this is side 3
					nodes_to_use = node_side3;
				}
				for (int j = 0; j < degree; j++) {
					point1[0] = triangles[tri_index]->controlP[nodes_to_use[j]][0];
					point1[1] = triangles[tri_index]->controlP[nodes_to_use[j]][1];

					point2[0] = triangles[tri_index]->controlP[nodes_to_use[j + 1]][0];
					point2[1] = triangles[tri_index]->controlP[nodes_to_use[j + 1]][1];

					section_len += sqrt(pow(point2[0] - point1[0], 2) + pow(point2[1] - point1[1], 2));   // update the distance
				}
				element_len += section_len;
				all_section_lens.push_back(section_len);
			}
			// Because I use the linear distance between triangular control points to determine the position within the bezier element,
			// I will use the sum of these distance instead of the curve lenght which was calculated in Part1 of this program.
			// These numbers would ideally be the same, however, their slight differences cause me to make this decision.

			// I will also take care of Step 3 within this loop (projecting the control points
			double last_xi = 0.0;

			for (unsigned int section = 0; section < same_elem_handle.size(); section++) {
				int tri_index = tri_NURB_elem_section_side[same_elem_handle[section]][0];

				// determine node indexing info (for step 3)
				vector<int> nodes_to_use;
				if (tri_NURB_elem_section_side[same_elem_handle[section]][4] == 0) {   // this is side 1
					nodes_to_use = node_side1;
				}
				else if (tri_NURB_elem_section_side[same_elem_handle[section]][4] == 1) {    // this is side 2
					nodes_to_use = node_side2;
				}
				else {     // this is side 3
					nodes_to_use = node_side3;
				}


				double section_portion = all_section_lens[section] / element_len;
				for (int j = 1; j <= degree + 1; j++) {
					double xi = ((double (j) / double (degree + 1)) * section_portion) + last_xi;
					xi_list.push_back(xi);

					// Project the control points.  I will only alter the point within triangles->controlP, not node_list
					// I will alter node_list only after all of the weight calculation is completed
					if (triangles[tri_index]->controlP[nodes_to_use[j - 1]][2] && triangles[tri_index]->controlP[nodes_to_use[j - 1]][2] != 1.0) {
						// do this if the weight is non zero.  This is because nonzero means I know what the weight is.  Also, skip the loop if the weight is exactly 1 since I don't really need to multiply by 1
						triangles[tri_index]->controlP[nodes_to_use[j - 1]][0] *= triangles[tri_index]->controlP[nodes_to_use[j - 1]][2];
						triangles[tri_index]->controlP[nodes_to_use[j - 1]][1] *= triangles[tri_index]->controlP[nodes_to_use[j - 1]][2];
					}
				}
				last_xi = xi_list.back();
			}
			/////////////////////////////////////////////////////////////
			/////////////////////   STEP 4   ////////////////////////////
			/////////////////////////////////////////////////////////////

			// Now I need to perform curve refinement and degree elevation.  Since this will be kind of lengthy, I will do it within another function

		}
	}
}

void nurb::adjust_boundary_deg3()
{
	{
		for (unsigned int i = 0; i < triangles.size(); i++) {  //  loop through all of the triangles
			for (int side = 0; side < 3; side++) {
				if (global_edges[triangles[i]->global_side[side]][4] != phy_groups) {  // if it passes this, it is a boundary edge
																								// now I need to loop through all of the nodes that would be on the boundary, there are degree + 1 of them
					bool found0 = false;
					bool found1 = false;
					bool found2 = false;
					bool found3 = false;  // these are for each of the 4 nodes along the edge

					for (unsigned int j = 0; j < bNodes.size(); j++) {
						if (bNodes[j] == global_edges[triangles[i]->global_side[side]][0]) {
							found0 = true;
						}
						if (bNodes[j] == global_edges[triangles[i]->global_side[side]][1]) {
							found1 = true;
						}
						if (bNodes[j] == global_edges[triangles[i]->global_side[side]][2]) {
							found2 = true;
						}
						if (bNodes[j] == global_edges[triangles[i]->global_side[side]][3]) {
							found3 = true;
						}
						if (found0 && found1 && found2 && found3)
							break;
					}
					if (found0 == false)  // since it was not found, add the node to bNodes
						bNodes.push_back(global_edges[triangles[i]->global_side[side]][0]);

					if (found1 == false)
						bNodes.push_back(global_edges[triangles[i]->global_side[side]][1]);

					if (found2 == false)
						bNodes.push_back(global_edges[triangles[i]->global_side[side]][2]);

					if (found3 == false)
						bNodes.push_back(global_edges[triangles[i]->global_side[side]][3]);
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////

					// now determine which of the KV sections ("half elements") this belongs to.  Currently, this is which half element this belongs to
					unsigned int KV_section = global_edges[triangles[i]->global_side[side]][4];   // KV_section starts at 1
					int cur_nurb = 0;   // this will be the identifier for the current nurbs curve.
					bool found = false;
					while (found == false) {
						// in the way the mesh is defined, there is a straight line connecting all of the xi_real_loc point for the NURBS curve

						unsigned xi_size = Elem_list[cur_nurb]->xi_evals.size();
						if (KV_section <= Elem_list[cur_nurb]->n_el * (xi_size - 1)) { // this means the KV_section lies in the current NURBS curve
							unsigned int element = KV_section / xi_size; // this is the identifier within the bezier element for that NURBS curve. Starts at 0
							unsigned int elem_section = KV_section % xi_size;    // this is the identifier for the section within that bezier element.  Starts at 0
							if (elem_section == 0) {
								element++;
								elem_section = xi_size - 2;     // the reason there is this -2 is because there is 1 fewer lines than there are xi points and need to subtact 1 for proper indexing
							}
							else
								elem_section--;
							// now that I know that stuff, I need determine what the parametric xi is which corresponds to the physical locations along the edge
							// first I need to determine the length of the line between the beginning and end of the "half element"
							vector<double> first_point(2, 0);
							vector<double> second_point(2, 0);
							first_point[0] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section][0];    // this is the x location for the first point
							first_point[1] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section][1];    // this is the y location for the first point

							second_point[0] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section + 1][0];    // this is the x location for the second point
							second_point[1] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section + 1][1];    // this is the y location for the second point

							double len = sqrt(pow(second_point[0] - first_point[0], 2) + pow(second_point[1] - first_point[1], 2));    // this is the linear length between the beginning and end of the "half element"

							vector <double> temp(2, 0); // this is a variable I will use to write the value of the current control point location into
							double temp_len = 0;
							for (int k = 0; k < 4; k++) {     //////////////////////////////////////////////////////////  I MIGHT WANT TO CHANGE THIS LOOP SO IT ONLY DOES THE CORNER NODE (K = 0,4)
								int node_index = global_edges[triangles[i]->global_side[side]][k];   // this will let me reference the correct node in the node list
								temp[0] = node_list[node_index][0];
								temp[1] = node_list[node_index][1];
								temp_len = sqrt(pow(temp[0] - first_point[0], 2) + pow(temp[1] - first_point[1], 2));
								temp_len = temp_len / len;
								// now I need to add in the parametric value of the beginning of this section.
								double begin_xi = Elem_list[cur_nurb]->xi_evals[elem_section];
								double end_xi = Elem_list[cur_nurb]->xi_evals[elem_section + 1];
								temp_len = ((end_xi - begin_xi) * temp_len) + begin_xi;   //  this is the value I need to evaluate the original bezier element at
								temp = eval_Bez_elem(temp_len, element, cur_nurb);
								// now I need to update the position of the node with in node_list
								node_list[node_index] = temp;
							}

							found = true;
						}
						else {
							KV_section -= Elem_list[cur_nurb]->n_el * (xi_size - 1);
							cur_nurb++;
						}
					}
				}
			}

		}
	}




}

void nurb::smooth_weights(int degree)
{


}





////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////  Deal with all of this stuff later//////////////////////
vector<double> nurb::eval_Bez_elem(double xi_val, unsigned int element, unsigned int cur_nurb)
{
	vector <double> num(2, 0);
	double denom = 0;
	Bezier_handle *Bez = Elem_list[cur_nurb];
	for (int k = 0; k < Bez->p + 1; k++) {            // Loop through all of the control points in that element.  There are p + 1 of them
													  // first I need to evaluate the basis functions at the temp_len
		double N = (fast_fact[Bez->p] / (fast_fact[k] * fast_fact[Bez->p - k]))*pow(xi_val, (k)) * pow((1 - xi_val), (Bez->p - k));
		num[0] += N * Bez->elem_Geom[element].controlP[k][0] * Bez->elem_Geom[element].weight[k];
		num[1] += N * Bez->elem_Geom[element].controlP[k][1] * Bez->elem_Geom[element].weight[k];
		denom += N * Bez->elem_Geom[element].weight[k];
	}
	num[0] = num[0] / denom;
	num[1] = num[1] / denom;


	return num;
}

void nurb::LE2D(int degree)
{
	// define the stiffness properties
	double E = pow(10, 11);
	double nu = 0.3;
	double lambda = (nu*E) / ((1 + nu)*(1 - 2 * nu));
	double mu = E / (2 * (1 + nu));
	Matrix3d D;
	D << lambda + 2 * mu, lambda, 0,
		lambda, lambda + 2 * mu, 0,
		0, 0, mu;

	// Mesh / Geometry Information
	const unsigned int nel = triangles.size();					// This is the number of triangular elements
	const unsigned int nen = 10;									// This is the number of nodes within each element
	const unsigned int nNodes = node_list.size();					// This is the number of unique nodes throughout the mesh
	const unsigned int nsd = 2;                                   // This is the number of spatial dimensions
	const unsigned int ndof = nsd*(nNodes - bNodes.size());       // This is the number of degrees of freedom
	const unsigned int nedof = nen * nsd;							// This is the number of element-wise degrees of freedom

																	// Initialize the global K and F matrices
	MatrixXd K(ndof, ndof);
	K.setZero();
	VectorXd F(ndof);
	F.setZero();

	// Now construct the ID array
	unsigned int P = 0;                  ////////// LUKE HAS THIS AS A 1,   FIND OUT WHICH ONE LATER
	MatrixXi ID(nsd, nNodes);
	ID.setZero();

	for (unsigned int a = 0; a < nNodes; a++) {
		////////// find if a exists with in bNodes   ///////////
		bool found = false;
		for (unsigned int j = 0; j < bNodes.size(); j++) {
			if (a == bNodes[j]) {
				found = true;
				break;
			}
		}
		////////////////////////////////////////////////////////

		if (found) {
			ID(0, a) = 0;
			ID(1, a) = 0;
		}
		else {
			ID(0, a) = P;
			ID(1, a) = P + 1;
			P += 2;
		}
	}

	// Generate a Lookup table for the basis
	evaluate_tri_basis(degree);

	// now loop through all of the elements
	for (unsigned int elem = 0; elem < triangles.size(); elem++) {
		// set the local stiffness matrix to 0
		MatrixXd k_loc(nedof, nedof);
		k_loc.setZero();

		// loop over quadrature points to construct k
		for (int q = 0; q < num_quad; q++) {
			// construct the B matrix
			MatrixXd B(3, nedof);
			B.setZero();
			tri_10_output fast = tri_10_fast(elem, q);
			cout << "block of size 5x2" << endl;
			for (unsigned int a = 0; a < nedof; a++) {

			}

		}

	}
}

void nurb::evaluate_tri_basis(int degree)
{
	int n = 3;
	quadInfo quad;
	// now I need to loop through all of these points

	// set up A matrix to so that barycentric coordinates can be solved for
	Matrix3d A;
	A << 0, 1, 0,
		0, 0, 1,
		1, 1, 1;
	Vector3d B;
	Vector3d bary;
	MatrixXi tuples(degree, 3);
	tuples << 3, 0, 0,
		0, 3, 0,
		0, 0, 3,
		2, 1, 0,
		1, 2, 0,
		0, 2, 1,
		0, 1, 2,
		1, 0, 2,
		2, 0, 1,
		1, 1, 1;

	for (int cur_point = 0; cur_point < num_quad; cur_point++) {    // This only goes up to 16 because there are only 16 quad points
																	// now determine the barycentric coordinate
		B << quad.quadP[cur_point][0], quad.quadP[cur_point][1], 1;
		bary = A.lu().solve(B);

		double u = bary(0);
		double v = bary(1);
		double w = bary(2);
		vector<double> temp(10, 0);   // this will hold each row of the tri_N before it is pushed back
		vector<vector<double>> temp_der;

		for (int m = 0; m < 10; m++) {
			int i = tuples(m, 0);
			int j = tuples(m, 1);
			int k = tuples(m, 2);

			temp[m] = fast_fact[n] / (fast_fact[i] * fast_fact[j] * fast_fact[k])*pow(u, i)*pow(v, j)*pow(w, k);
			vector <double> deriv(3, 0);
			if ((i - 1) >= 0)
				deriv[0] = n * fast_fact[n - 1] / (fast_fact[i - 1] * fast_fact[j] * fast_fact[k]) * pow(u, i - 1) * pow(v, j) * pow(w, k);
			if ((j - 1) >= 0)
				deriv[1] = n * fast_fact[n - 1] / (fast_fact[i] * fast_fact[j - 1] * fast_fact[k]) * pow(u, i) * pow(v, j - 1) * pow(w, k);
			if ((k - 1) >= 0)
				deriv[2] = n * fast_fact[n - 1] / (fast_fact[i] * fast_fact[j] * fast_fact[k - 1]) * pow(u, i) * pow(v, j) * pow(w, k - 1);

			temp_der.push_back(deriv);

		}
		tri_N.push_back(temp);
		tri_dN_du.push_back(temp_der);
	}
}
*/


//tri_10_output nurb::tri_10_fast(unsigned int elem, int q)
//{

	/*// testing these are the same nodes as for the first call to tri10fast with plateandhole
	triangles[elem]->controlP[0][0] = -4.00000;
	triangles[elem]->controlP[0][1] = -4.00000;
	triangles[elem]->controlP[1][0] = -2.00000;
	triangles[elem]->controlP[1][1] = -4.00000;
	triangles[elem]->controlP[2][0] = -2.6919;
	triangles[elem]->controlP[2][1] = -2.6919;
	triangles[elem]->controlP[3][0] = -3.3333;
	triangles[elem]->controlP[3][1] = -4.0000;
	triangles[elem]->controlP[4][0] = -2.6667;
	triangles[elem]->controlP[4][1] = -4.0000;
	triangles[elem]->controlP[5][0] = -2.2306;
	triangles[elem]->controlP[5][1] = -3.5640;
	triangles[elem]->controlP[6][0] = -2.4613;
	triangles[elem]->controlP[6][1] = -3.1279;
	triangles[elem]->controlP[7][0] = -3.1279;
	triangles[elem]->controlP[7][1] = -3.1279;
	triangles[elem]->controlP[8][0] = -3.5640;
	triangles[elem]->controlP[8][1] = -3.5640;
	triangles[elem]->controlP[9][0] = -2.8973;
	triangles[elem]->controlP[9][1] = -3.5640;

	*/
	/*

		unsigned int nen = tri_N[0].size();
		tri_10_output output;   // create instance of struct I will eventually return
		output.R.resize(nen);

		MatrixXd dR_du(nen, 3);
		dR_du.setZero();
		// to save on loop unrolling I will hard code all of this out.
		// the next three variables are the sums of tri_dN_du along each of the second dimensions.
		// this is equivilent to
		// dN_du(:,1)'*node(:,3) in line 59 of tri10fast of Luke's Matlab code
		double dN_du0_sum = tri_dN_du[q][0][0]
			+ tri_dN_du[q][1][0]
			+ tri_dN_du[q][2][0]
			+ tri_dN_du[q][3][0]
			+ tri_dN_du[q][4][0]
			+ tri_dN_du[q][5][0]
			+ tri_dN_du[q][6][0]
			+ tri_dN_du[q][7][0]
			+ tri_dN_du[q][8][0]
			+ tri_dN_du[q][9][0];

		double dN_du1_sum = tri_dN_du[q][0][1]
			+ tri_dN_du[q][1][1]
			+ tri_dN_du[q][2][1]
			+ tri_dN_du[q][3][1]
			+ tri_dN_du[q][4][1]
			+ tri_dN_du[q][5][1]
			+ tri_dN_du[q][6][1]
			+ tri_dN_du[q][7][1]
			+ tri_dN_du[q][8][1]
			+ tri_dN_du[q][9][1];

		double dN_du2_sum = tri_dN_du[q][0][2]
			+ tri_dN_du[q][1][2]
			+ tri_dN_du[q][2][2]
			+ tri_dN_du[q][3][2]
			+ tri_dN_du[q][4][2]
			+ tri_dN_du[q][5][2]
			+ tri_dN_du[q][6][2]
			+ tri_dN_du[q][7][2]
			+ tri_dN_du[q][8][2]
			+ tri_dN_du[q][9][2];

		vector<double>num(2, 0.0);
		double den = 0.0;

		for (unsigned int i = 0; i < nen; i++) {
			output.R[i] = tri_N[q][i];

			dR_du(i, 0) = tri_dN_du[q][i][0] - (dN_du0_sum * tri_N[q][i]);
			dR_du(i, 1) = tri_dN_du[q][i][1] - (dN_du1_sum * tri_N[q][i]);
			dR_du(i, 2) = tri_dN_du[q][i][2] - (dN_du2_sum * tri_N[q][i]);

			num[0] += tri_N[q][i] * triangles[elem]->controlP[i][0];
			num[1] += tri_N[q][i] * triangles[elem]->controlP[i][1];
			den += tri_N[q][i];

		}
		num[0] = num[0] / den;
		num[1] = num[1] / den;
		output.x = num;

		//////////////////////////////////////////////////////////
		// now Chain rule to find the derivative with respect to cartesian isoparametric coordinates
		MatrixXd du_dxi(3, 2);
		du_dxi << -1, -1,
			1, 0,
			0, 1;

		MatrixXd dR_dxi = dR_du * du_dxi;
		MatrixXd dN_du_mat = create_matrix(tri_dN_du[q]);
		MatrixXd dN_dxi = dN_du_mat * du_dxi;

		// now calculate the mapping from isoparametric space to physical space
		Matrix2d g, gp, h, hp;
		h.setConstant(1);
		g.setZero();
		gp.setZero();
		hp.setZero();
		for (int row = 0; row < 2; row++) {
			for (int col = 0; col < 2; col++) {
				for (unsigned int i = 0; i < nen; i++) {
					g(row, col) += tri_N[q][i] * triangles[elem]->controlP[i][row];
					gp(row, col) += dN_dxi(i, col) * triangles[elem]->controlP[i][row];
					hp(row, col) += dN_dxi(i, col);
				}
			}
		}
		output.Jacob = (gp.array() * h.array()) - (g.array() * hp.array());      // determine the jacobian of the mapping
		output.Jacob = output.Jacob.matrix();

		MatrixXd temp(2, 2);
		temp << output.Jacob.inverse();
		output.dR_dx = dR_dxi * temp;
		output.J_det = output.Jacob.determinant();

		cout << output.dR_dx << endl << endl;
		cout << output.Jacob << endl << endl;
		cout << output.J_det << endl << endl;


		return output;
	}

	MatrixXd nurb::create_matrix(std::vector<std::vector<double>> input)
	{
		unsigned first = input.size();
		unsigned second = input[0].size();
		MatrixXd output(first, second);

		for (unsigned int i = 0; i < first; i++) {
			for (unsigned int j = 0; j < second; j++) {
				output(i, j) = input[i][j];
			}
		}
		return output;
	}
	*/
