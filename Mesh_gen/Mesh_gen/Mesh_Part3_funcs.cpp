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
#include <Eigen/SparseLU>
#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

void nurb::smoothMesh(int mesh_degree)
{
	organize_boundary();
	split_and_extract();

	// now I need to smooth the weights of all of the points
	smooth_weights(degree);


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

					// in the way the mesh is defined, there is a straight line connecting all of the xi_real_loc point for the NURBS curve

					unsigned int xi_size = int(Elem_list[cur_nurb]->xi_evals.size());
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
	size_t size = tri_NURB_elem_section_side.size();
	vector<int> work_around(size);
	iota(begin(work_around), end(work_around), 0);
	sort(work_around.begin(), work_around.end(), *this);
	vector<vector<unsigned int>> temp(size);
	for (unsigned int i = 0; i < size; i++) {
		temp[i] = tri_NURB_elem_section_side[work_around[i]];
	}
	tri_NURB_elem_section_side = temp;

	// now go through and add an identifier as to whether or not the side is "backwards"
	int count = 0;
	for (int i = 0; i < num_curves; i++) {
		for (int j = 0; j < Elem_list[i]->n_el; j++) {
			vector<double> anchor = Elem_list[i]->elem_Geom[j].controlP[0];

			int elem = tri_NURB_elem_section_side[count][2];

			while (count < int(tri_NURB_elem_section_side.size()) && tri_NURB_elem_section_side[count][2] == elem) {

				int tri = tri_NURB_elem_section_side[count][0];
				int side = tri_NURB_elem_section_side[count][4];
				int node = node_side_index[side][0];

				vector<double> proper(2, 0);

				int node_list_index = triangles[tri]->controlP[node];
				proper[0] = node_list[node_list_index][0];
				proper[1] = node_list[node_list_index][1];
				double dist_prop = sqrt(pow(proper[0] - anchor[0], 2) + pow(proper[1] - anchor[1], 2));

				vector<double> reverse(2, 0);
				int rev_node = node_side_index[side][degree];
				int rev_node_list_index = triangles[tri]->controlP[rev_node];

				reverse[0] = node_list[rev_node_list_index][0];
				reverse[1] = node_list[rev_node_list_index][1];
				double dist_rev = sqrt(pow(reverse[0] - anchor[0], 2) + pow(reverse[1] - anchor[1], 2));

				if (dist_prop < dist_rev) {  // this means that it is in the "correct" order
					tri_NURB_elem_section_side[count].push_back(1);
					anchor = reverse;
				}
				else {   // this means that it in the reverse order
					tri_NURB_elem_section_side[count].push_back(0);
					anchor = proper;
				}

				count++;
			}
		}
	}
	// now organize bNodes to that it is an increasing list
	sort(bNodes.begin(), bNodes.end());
	// now loop through to create bNodeBool
	bNodeBool.resize(node_list.size());
	count = 0;
	for (int i = 0; i < int(node_list.size()); i++) {
		if (bNodes[count] == i) {
			bNodeBool[i] = true;
			count++;
			if (count == bNodes.size())
				break;
		}
	}
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
		else {   // this means that the element indexes are the same
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
			int num_to_increment = int(xi_to_add.size());
			// now I need to add enough knots at these locations so that there is c(0) continuity
			int nurb_degree = Elem_list[cur_nurb]->p;
			for (unsigned int i = 0; i < xi_to_add.size(); i++) {
				vector<double> repeats(nurb_degree - 1, xi_to_add[i]);
				xi_to_add.insert(xi_to_add.begin() + i, repeats.begin(), repeats.end());
				i += int(repeats.size());
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
	vector<double> split_list;

	unsigned int max_i = int(tri_NURB_elem_section_side.size());
	while (unsigned(i) < max_i && tri_NURB_elem_section_side[i][2] == cur_elem && tri_NURB_elem_section_side[i][1] == cur_nurb) {
		// do Newton's method to find the exact parametric xi to split at
		// need to figure out a way to make sure that second_point is actually the second point along the NURBS curve
		vector<double> second_point;
		if (tri_NURB_elem_section_side[i][5])  // the edge is in the proper direction
			second_point = node_list[triangles[tri_NURB_elem_section_side[i][0]]->controlP[node_side_index[tri_NURB_elem_section_side[i][4]][degree]]];
		else   // the edge is in the reverse direction
			second_point = node_list[triangles[tri_NURB_elem_section_side[i][0]]->controlP[node_side_index[tri_NURB_elem_section_side[i][4]][0]]];

		double value = newton_find_xi(cur_nurb, cur_elem, second_point, 0.5);
		split_list.push_back(value);
		i++;
	}
	split_list.pop_back();
	return split_list;
}

void nurb::curve_refine(Bezier_handle * Bez, int cur_elem, vector<double> xi_to_add, bool part1)
{
	// Step 1: create a knot vector for this element.  Since it is Bezier, it is easy to compute
	int p = Bez->p;
	unsigned int n = int(Bez->elem_Geom[cur_elem].controlP.size());
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
		unsigned int k = int(index - 1);
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
		if (i) {
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
				if (cur_elem + j >= Bez->Operator.size())
					Bez->Operator.push_back(Bez->Operator[cur_elem + j - 1]);
				else
					Bez->Operator.insert(Bez->Operator.begin() + cur_elem + j, Bez->Operator[cur_elem + j]);

				// Now Update the master IEN array so it reflects the added control points.
				vector<int> insert(Bez->p + 1, 0);  // this is the entry corresponding to the inserted point
				vector<vector<int>> cur_IEN = Master_IEN.back();
				insert[0] = cur_IEN[cur_elem + j - 1][Bez->p];
				for (int i = 1; i < Bez->p + 1; i++) {
					insert[i] = insert[i - 1] + 1;
				}
				if (cur_elem + j >= cur_IEN.size())
					cur_IEN.push_back(insert);
				else
					cur_IEN.insert(cur_IEN.begin() + cur_elem + j, insert);

				for (int k = cur_elem + 1 + j; k < int(cur_IEN.size()); k++) {   
					for (int m = 0; m < Bez->p + 1; m++) {
						cur_IEN[k][m] += Bez->p;
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
	int p = int(Bez->elem_Geom[element].controlP.size()) - 1;
	unsigned int n = int(P.size());
	Pe.push_back(Pe[n - 1]);

	for (unsigned int i = 1; i < n; i++) {
		Pe[i][0] = (double(i) / double(p + 1))*P[i - 1][0] + (1 - (double(i) / double(p + 1)))*P[i][0];
		Pe[i][1] = (double(i) / double(p + 1))*P[i - 1][1] + (1 - (double(i) / double(p + 1)))*P[i][1];
		Pe[i][2] = (double(i) / double(p + 1))*P[i - 1][2] + (1 - (double(i) / double(p + 1)))*P[i][2];

	}

	Bez->elem_Geom[element].controlP = Pe;

	// now figure out what to make the degree
	unsigned int min_p = int(Bez->elem_Geom[0].controlP.size()) - 1;
	for (int i = 1; i < Bez->n_el; i++) {
		if (Bez->elem_Geom[i].controlP.size() - 1 < min_p)
			min_p = int(Bez->elem_Geom[i].controlP.size()) - 1;
	}
	Bez->p = min_p;

}

vector<vector<double>> nurb::determine_dirichlet()
{
	vector<vector<double>> g(node_list.size(),vector<double>(3,0));
	// loop through tri_NURB_elem_section_side
	for (unsigned int i = 0; i < tri_NURB_elem_section_side.size(); i++) {

		int tri = tri_NURB_elem_section_side[i][0];
		int nurb = tri_NURB_elem_section_side[i][1];
		int elem = tri_NURB_elem_section_side[i][2];
		int side = tri_NURB_elem_section_side[i][4];
		int order = tri_NURB_elem_section_side[i][5];

		// loop through all of the points except for the first point'
		// since this list is in order and these curves will always form a closed loop, I can just do the last p points instead of all p+1 of them
		for (int j = 1; j <= degree; j++) {
			// determine the correct value of the control point
			int node = node_side_index[side][j];

			vector<double> actual_p(3, 0);
			actual_p[0] = Elem_list[nurb]->elem_Geom[elem].controlP[j][0];
			actual_p[1] = Elem_list[nurb]->elem_Geom[elem].controlP[j][1];
			actual_p[2] = Elem_list[nurb]->elem_Geom[elem].controlP[j][2];


			vector<double> g_row(3, 0);
			// fill g_row with the displacement dirichlet condition
			// At this point, all of the the weights in the node_list will be 1

			if (order) { // this is in the proper direction
				// determine which node index this is
				int node_list_index = triangles[tri]->controlP[node];
				g_row[0] = actual_p[0] - node_list[node_list_index][0];
				g_row[1] = actual_p[1] - node_list[node_list_index][1];
				g_row[2] = actual_p[2] - node_list[node_list_index][2];
				g[node_list_index] = g_row;
			}
			else {   // this is in the reverse direction
				int rev_node = node_side_index[side][degree - j];
				int rev_node_list_index = triangles[tri]->controlP[rev_node];
				g_row[0] = actual_p[0] - node_list[rev_node_list_index][0];
				g_row[1] = actual_p[1] - node_list[rev_node_list_index][1];
				g_row[2] = actual_p[2] - node_list[rev_node_list_index][2];
				g[rev_node_list_index] = g_row;
			}
		}
	}
	// now determine the neighbor_LE_nodes to help speed up the element formation
	int list_size = (int(node_list.size()) - int(bNodes.size())) << 1;
	neighbor_LE_num.resize(list_size);
	int count = 0;  // this will keep track of its position in the unaltered list
	for (int i = 0; i < list_size; i++) {
		// determine if the entire row needs to be erased or not
		if (bNodeBool[count]) {  // the row needs to be deleted
			neighbor_nodes.erase(neighbor_nodes.begin() + i);
			i--;
		}
		else { // we need to go down the length of the row and determine if the any nodes need to be deleted from within
			int row_len = int(neighbor_nodes[i].size());
			for (int j = 0; j < row_len; j++) {
				if (bNodeBool[neighbor_nodes[i][j]]) {   // it needs to be deleted
					neighbor_nodes[i].erase(neighbor_nodes[i].begin() + j);
					j--;
					row_len--;
				}
				else {  // it needs to be altered
					neighbor_nodes[i][j] = neighbor_nodes[i][j] << 1;
					neighbor_nodes[i].insert(neighbor_nodes[i].begin() + j + 1, neighbor_nodes[i][j] + 1);
					j++;
					row_len++;
				}
			}
			// I need to make a duplicate row
			neighbor_nodes.insert(neighbor_nodes.begin() + i + 1, neighbor_nodes[i]);
			neighbor_LE_num(i) = row_len;
			neighbor_LE_num(i + 1) = row_len;
			i++;
		}
		count++;
	}
	return g;
}

double nurb::newton_find_xi(int cur_nurb, int cur_elem, vector<double> second_point, double guess)
{
	double tol = 0.0000001;
	double error = 1.0;
	int p = Elem_list[cur_nurb]->p;
	int iter = 0;     // set a max number of iteration to prevent infinite loops
	while (error > tol) {   // carry out the newton's method
		// find the slope of the line, this is just the derivative of the bezier coordinate since this is the only part which changes with time
		if (iter > 15)
			return guess;


		vector<double> slope_vec(2, 0);
		vector<double> act_loc(2, 0);
		//double bez_der_tot = 0;
		double bez_tot = 0;
		vector<double> bez_part_der(p + 1, 0);
		vector<double> bez_part(p + 1, 0);

		for (int i = 0; i <= p; i++) {   // loop through all of the control point throughout the bezier element
			vector<double> point = Elem_list[cur_nurb]->elem_Geom[cur_elem].controlP[i];

			bez_part_der[i] = ((n_choose_k(p, i) * i * pow(guess, i - 1) * pow((1 - guess), (p - i))) + (n_choose_k(p, i) * pow(guess, i) * (i - p) * pow((1 - guess), (p - i - 1))));
			//bez_der_tot = bez_der_tot + (bez_part_der[i] * point[2]);

			bez_part[i] = n_choose_k(p, i) * pow(guess, i) * pow((1 - guess), (p - i));
			bez_tot = bez_tot + (bez_part[i] * point[2]);
		}

		for (int i = 0; i <= p; i++) {
			vector<double> point = Elem_list[cur_nurb]->elem_Geom[cur_elem].controlP[i];

			// deal with the weights
			//bez_part_der[i] = (bez_part_der[i] * point[2]) / bez_der_tot;
			bez_part[i] = (bez_part[i] * point[2]) / bez_tot;

			slope_vec[0] = slope_vec[0] + point[0] * bez_part_der[i];
			slope_vec[1] = slope_vec[1] + point[1] * bez_part_der[i];
			act_loc[0] = act_loc[0] + point[0] * bez_part[i];
			act_loc[1] = act_loc[1] + point[1] * bez_part[i];
		}


		// now determine what x(xi) - x_tilde is
		error = sqrt(pow(act_loc[0] - second_point[0], 2) + pow(act_loc[1] - second_point[1], 2));  // this is just the distance formula
		double slope = (((act_loc[0] - second_point[0]) * slope_vec[0]) + ((act_loc[1] - second_point[1]) * slope_vec[1])) / error;

		// return what guess is now if there is no error
		if (!error)
			return guess;

		// now determine where this line crosses the zero
		guess = guess - (error / slope);
		iter++;
	}
	return guess;
}

double nurb::n_choose_k(int n, int k)
{
	double ans;
	ans = fast_fact[n] / (fast_fact[k] * fast_fact[n - k]);
	return ans;
}

void nurb::smooth_weights(int degree)
{

	vector<vector<double>> g = determine_dirichlet();

	// first I need to smooth the weights by solving the laplacian
	solve_laplacian(g);

	// now solve the linear elasticity problem
	LE2D(g);


}

void nurb::solve_laplacian(vector<vector<double>> g)
{
	// Declare the Global K and F matrices
	int nNodes = int(node_list.size());
	SparseMatrix<double, ColMajor> K(nNodes, nNodes);
	K.reserve(neighbor_nodes_num);
	VectorXd F(nNodes);
	F.setZero();

	evaluate_tri_basis(28);
	quadInfo28 info;

	MatrixXd k_loc(nodes_in_triangle, nodes_in_triangle);
	for (unsigned int i = 0; i < triangles.size(); i++) {
		// set the local stiffness matrix to 0
		k_loc.setZero();

		for (unsigned int q = 0; q < tri_N.size(); q++) {
			tri_10_output fast = tri_10_fast(i, q, false);

			fast.J_det = abs(fast.J_det);
			k_loc = k_loc + info.weights[q] * fast.dR_dx * fast.dR_dx.transpose() * fast.J_det;

		}
		// now assemble the global K and F matrices
		for (int b = 0; b < nodes_in_triangle; b++) {
			for (int a = 0; a <= b; a++) {
				K.coeffRef(triangles[i]->controlP[a], triangles[i]->controlP[b]) += k_loc(a, b);
			}
		}
	}
	// now fill the the bottom part of the K matrix
	SparseMatrix<double, ColMajor> K_trans = K.transpose();
	SparseMatrix<double, ColMajor> Diag = K;
	struct keep_diag {
		inline bool operator() (const int& row, const int& col, const double&) const
		{
			return row == col;
		}
	};
	Diag.prune(keep_diag());
	K += K_trans - Diag;

	// now take care of the boundary conditions
	// first match up each member of bNodes with the correct NURBS curve point
	vector<double> correct_weights(bNodes.size());
	for (unsigned int i = 0; i < bNodes.size(); i++) {
		int curve = 0;
		int elem = 0;
		int point = 0;
		int j = 0;
		bool found = false;
		while (!found) {  // it is gaurenteed to find it so I can do this and hit a break later
			Tri_elem *current = triangles[tri_NURB_elem_section_side[j][0]];
			curve = tri_NURB_elem_section_side[j][1];
			elem = tri_NURB_elem_section_side[j][2];
			int side = tri_NURB_elem_section_side[j][4];
			// the section doesn't matter since there is only one section per element at this point

			for (int k = 0; k < 3 * degree; k++) {
				if (current->controlP[k] == bNodes[i]) { // then we have a match
					// now use node_side_index to determine the node along the nurbs curve
					for (int m = 0; m <= degree; m++) {
						if (node_side_index[side][m] == k) {
							point = m;
							found = true;
							break;
						}
					}
					break;
				}
			}
			// increment the j for the while loop
			j++;
		}
		correct_weights[i] = Elem_list[curve]->elem_Geom[elem].controlP[point][2];
	}

	// now update the F vector so that it has the correct values
	for (unsigned int i = 0; i < bNodes.size(); i++) {

		int col = bNodes[i];
		SparseMatrix<double, ColMajor> cur_col = K.col(col);
		int *non_zero_rows = cur_col.innerIndexPtr();
		for (int j = 0; j < cur_col.nonZeros(); j++) {
			int row = non_zero_rows[j];
			if (row != col)
				F(row) = F(row) - K.coeff(row, col) * correct_weights[i];
		}
		// now I need to make zeros in the other direction
		SparseMatrix<double, ColMajor> cur_row = K.row(col);

		for (int j = 0; j < cur_col.nonZeros(); j++) {
			int row = non_zero_rows[j];
			K.coeffRef(row, col) = 0.0;
			K.coeffRef(col, row) = 0.0;
		}

		// now set the ones along the diagonal
		K.coeffRef(col, col) = 1;
		F(col) = correct_weights[i];

	}

	K.prune(0.0);

	// now solve the system
	SparseLU<SparseMatrix<double>> solver;
	solver.analyzePattern(K);
	solver.factorize(K);

	VectorXd d = solver.solve(F);   // this is a SparseLU solver
	// now update the node_list with this new smoothed weights
	for (int i = 0; i < nNodes; i++) {
		node_list[i][2] = d(i);
	}
}

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

void nurb::LE2D(vector<vector<double>> g)
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
	const unsigned int nel = int(triangles.size());					// This is the number of triangular elements
	const unsigned int nen = nodes_in_triangle;									// This is the number of nodes within each element
	const unsigned int nNodes = int(node_list.size());					// This is the number of unique nodes throughout the mesh
	const unsigned int nsd = 2;                                   // This is the number of spatial dimensions                          THIS WOULD NEED TO BE 3 FOR PROJECTED
	const unsigned int ndof = nsd*(nNodes - int(bNodes.size()));       // This is the number of degrees of freedom
	const unsigned int nedof = nen * nsd;							// This is the number of element-wise degrees of freedom

																	// Initialize the global K and F matrices
	SparseMatrix<double, ColMajor> K(ndof, ndof);
	K.reserve(neighbor_LE_num);
	VectorXd F(ndof);
	F.setZero();

	// Now construct the ID array
	unsigned int P = 0;                  ////////// LUKE HAS THIS AS A 1,   FIND OUT WHICH ONE LATER
	vector<vector<int>> ID(nNodes, vector<int>(nsd, -1));

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
			ID[a][0] = -1;
			ID[a][1] = -1;
		}
		else {
			ID[a][0] = P;
			ID[a][1] = P + 1;
			P += 2;
		}
	}
	// Generate a Lookup table for the basis

	tri_N.erase(tri_N.begin(), tri_N.end());
	tri_dN_du.erase(tri_dN_du.begin(), tri_dN_du.end());
	evaluate_tri_basis(16);
	quadInfo info;
	// now loop through all of the elements
	for (unsigned int elem = 0; elem < triangles.size(); elem++) {

		// set the local stiffness matrix to 0
		MatrixXd k_loc(nedof, nedof);
		k_loc.setZero();

		// loop over quadrature points to construct k
		for (int q = 0; q < num_quad; q++) {
			// construct the B matrix
			MatrixXd B(3, nedof);
			tri_10_output fast = tri_10_fast(elem, q, true);
			for (unsigned int a = 0; a < nen; a++) {
				MatrixXd insert_mat(3, 2);
				insert_mat << fast.dR_dx(a, 0), 0, 0, fast.dR_dx(a, 1), fast.dR_dx(a, 1), fast.dR_dx(a, 0);
				B.block(0, nsd*a, 3, 2) = insert_mat;
			}
			fast.J_det = abs(fast.J_det);
			k_loc = k_loc + info.weights[q]*B.transpose()*D*B*fast.J_det;
		}
		
		// assemble the local element stiffness matrix to the global stiffness matrix
		for (unsigned int a = 0; a < nen; a++) {
			for (unsigned int b = 0; b < nen; b++) {
				for (unsigned int i = 0; i < nsd; i++) {
					for (unsigned int j = 0; j < nsd; j++) {
						unsigned int p = a*nsd + i;
						unsigned int q = b*nsd + j;
						unsigned int Pa = ID[triangles[elem]->controlP[a]][i];
						unsigned int Pb = ID[triangles[elem]->controlP[b]][j];

						if (Pa != (-1) && Pb != (-1)) {
							K.coeffRef(Pa, Pb) = K.coeff(Pa, Pb) + k_loc(p, q);
						}
						else if (Pa != (-1) && Pb == (-1)) {
							F(Pa) = F(Pa) - (k_loc(p, q) * g[triangles[elem]->controlP[b]][j]);
						}
					}
				}
			}
		}
	}


	// Solve the system
	SparseLU<SparseMatrix<double>> solver;
	K.makeCompressed();
	solver.analyzePattern(K);
	solver.factorize(K);
	VectorXd d = solver.solve(F);   // this is a SparseLU solver
	// now update the node list
	for (unsigned int i = 0; i < nNodes; i++) {

		int xidx = ID[i][0];
		int yidx = ID[i][1];
		if (xidx != (-1)) {
			node_list[i][0] += d(xidx);
			node_list[i][1] += d(yidx);
		}
		node_list[i][0] += g[i][0];
		node_list[i][1] += g[i][1];
	}


	//// now normalize the minimum determinant against the area of the triangle
	//for (unsigned int elem = 0; elem < triangles.size(); elem++) {

	//	double elem_det = numeric_limits<double>::max();
	//	//double elem_det = 0;
	//	for (int q = 0; q < num_quad; q++) {
	//		// construct the B matrix
	//		tri_10_output fast = tri_10_fast(elem, q, true);
	//		if (abs(fast.J_det) < elem_det)
	//			elem_det = fast.J_det;
	//	}

	//	Matrix3d Area;
	//	Area << node_list[triangles[elem]->controlP[0]][0], node_list[triangles[elem]->controlP[1]][0], node_list[triangles[elem]->controlP[2]][0],
	//		node_list[triangles[elem]->controlP[0]][1], node_list[triangles[elem]->controlP[1]][1], node_list[triangles[elem]->controlP[2]][1],
	//		1, 1, 1;
	//	double A = abs(0.5 * Area.determinant());
	//	norm_dets.push_back(elem_det / A);
	//	cout << norm_dets.back() << endl;
	//}
}

void nurb::evaluate_tri_basis(int q_points)
{
	int n = degree;
	if (q_points == 16) {

		quadInfo quad;
		// now I need to loop through all of these points

		// set up A matrix to so that barycentric coordinates can be solved for
		Matrix3d A;
		A << 0, 1, 0,
			0, 0, 1,
			1, 1, 1;
		Vector3d B;
		Vector3d bary;

		for (int cur_point = 0; cur_point < num_quad; cur_point++) {    // This only goes up to 16 because there are only 16 quad points
																		// now determine the barycentric coordinate
			B << quad.quadP[cur_point][0], quad.quadP[cur_point][1], 1;
			bary = A.lu().solve(B);

			double u = bary(0);
			double v = bary(1);
			double w = bary(2);
			vector<double> temp(nodes_in_triangle, 0);   // this will hold each row of the tri_N before it is pushed back
			vector<vector<double>> temp_der;

			for (int m = 0; m < nodes_in_triangle; m++) {
				int i = int(bary_template[m][0] * double(degree));
				int j = int(bary_template[m][1] * double(degree));
				int k = int(bary_template[m][2] * double(degree));


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

	// NOW INCLUDE THE OPTION FOR THE 28 POINT QUADRATURE RULE
	else if(q_points == 28){
		quadInfo28 quad;
		// now I need to loop through all of these points

		// set up A matrix to so that barycentric coordinates can be solved for
		Matrix3d A;
		A << 0, 1, 0,
			0, 0, 1,
			1, 1, 1;
		Vector3d B;
		Vector3d bary;

		for (int cur_point = 0; cur_point < 28; cur_point++) {    // This only goes up to 16 because there are only 16 quad points
																		// now determine the barycentric coordinate
			B << quad.quadP[cur_point][0], quad.quadP[cur_point][1], 1;
			bary = A.lu().solve(B);

			double u = bary(0);
			double v = bary(1);
			double w = bary(2);
			vector<double> temp(nodes_in_triangle, 0);   // this will hold each row of the tri_N before it is pushed back
			vector<vector<double>> temp_der;

			for (int m = 0; m < nodes_in_triangle; m++) {
				int i = int(bary_template[m][0] * double(degree));
				int j = int(bary_template[m][1] * double(degree));
				int k = int(bary_template[m][2] * double(degree));


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
}

tri_10_output nurb::tri_10_fast(unsigned int tri, int q, bool rational)
{
	unsigned int nen = int(tri_N[0].size());
	tri_10_output output;   // create instance of struct I will eventually return
	output.R.resize(nen);

	MatrixXd dR_du(nen, 3);
	dR_du.setZero();
	

	// figure out and what all of the summed variables are
	// based on where I am in the program, I can ensure that "den" from Luke's code will always be 1
	// also, all of the weights will be 1 as well
	double dN_du0_sum = 0.0;
	double dN_du1_sum = 0.0;
	double dN_du2_sum = 0.0;


	for (unsigned int i = 0; i < nen; i++) {
		dN_du0_sum += tri_dN_du[q][i][0];
		dN_du1_sum += tri_dN_du[q][i][1];
		dN_du2_sum += tri_dN_du[q][i][2];

	}

	vector<double>num(3, 0.0);
	double den = 0;
	
	for (unsigned int i = 0; i < nen; i++) {
		if (rational) {

			dR_du(i, 0) = tri_dN_du[q][i][0] - (dN_du0_sum * tri_N[q][i]);
			dR_du(i, 1) = tri_dN_du[q][i][1] - (dN_du1_sum * tri_N[q][i]);
			dR_du(i, 2) = tri_dN_du[q][i][2] - (dN_du2_sum * tri_N[q][i]);
		}

		else {

			dR_du(i, 0) = tri_dN_du[q][i][0];
			dR_du(i, 1) = tri_dN_du[q][i][1];
			dR_du(i, 2) = tri_dN_du[q][i][2];
		}
		output.R[i] = tri_N[q][i];

		num[0] += tri_N[q][i] * node_list[triangles[tri]->controlP[i]][0];
		num[1] += tri_N[q][i] * node_list[triangles[tri]->controlP[i]][1];
		den += tri_N[q][i];

	}
	num[0] = num[0] / den;
	num[1] = num[1] / den;

	output.x = num;

	//////////////////////////////////////////////////////////
	// now Chain rule to find the derivative with respect to cartesian isoparametric coordinates
	MatrixXd du_dxi(3, 2);
	du_dxi << -1, -1,
		        1, 0,         // this might not be right if I need to go into projected space
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
				g(row, col) += tri_N[q][i] * node_list[triangles[tri]->controlP[i]][row];
				gp(row, col) += dN_dxi(i, col) * node_list[triangles[tri]->controlP[i]][row];
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

	return output;
}

MatrixXd nurb::create_matrix(std::vector<std::vector<double>> input)
{
	unsigned first = unsigned int(input.size());
	unsigned second = unsigned int(input[0].size());
	MatrixXd output(first, second);

	for (unsigned int i = 0; i < first; i++) {
		for (unsigned int j = 0; j < second; j++) {
			output(i, j) = input[i][j];
		}
	}
	return output;
}
