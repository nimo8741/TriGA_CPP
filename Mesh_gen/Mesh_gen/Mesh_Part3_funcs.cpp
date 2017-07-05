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

void nurb::smoothMesh(int degree)
{
	// first I need to find all of the boundary triangles and adjust all of the boundary control points so the line exactly on the geometry
	//if (degree == 63) {
	//	adjust_boundary_deg3();
	//}
	//else
		adjust_boundary(degree);

	boundary_weights(degree);
	// now I need to solve for the weights of all of the points]
	smooth_weights(degree);

	// now I need to solve the linear elasticity
	LE2D(degree);
}

void nurb::adjust_boundary(int degree)
{
	for (unsigned int i = 0; i < triangles.size(); i++) {  //  loop through all of the triangles
		for (int side = 0; side < 3; side++) {
			if (global_edges[triangles[i]->global_side[side]][degree + 1] != phy_groups) {  // if it passes this, it is a boundary edge
				// now I need to loop through all of the nodes that would be on the boundary, there are degree + 1 of them

				for (int edge_node = 0; edge_node < degree + 1; edge_node++) { // loop through all of the boundary nodes
					bool found = false;
					for (unsigned int j = 0; j < bNodes.size(); j++) {
						if (bNodes[j] == global_edges[triangles[i]->global_side[side]][edge_node]) {
							found = true;
							break;
						}
					}
					// if the node was not found within bNodes, add it to the list
					if (!found) {
						bNodes.push_back(global_edges[triangles[i]->global_side[side]][edge_node]);
					}


					// now determine which of the KV sections ("half elements") this belongs to.  Currently, this is which half element is the one denoted by the last number in global_edges
					unsigned int KV_section = global_edges[triangles[i]->global_side[side]][degree + 1];   // KV_section starts at 1
					unsigned int cur_nurb = 0;   // this will be the identifier for the current nurbs curve.
					found = false;
					while (found == false) {
						// in the way the mesh is defined, there is a straight line connecting all of the xi_real_loc point for the NURBS curve
													
						unsigned xi_size = Elem_list[cur_nurb]->xi_evals.size();
						if (KV_section < Elem_list[cur_nurb]->n_el * (xi_size - 1)) { // this means the KV_section lies in the current NURBS curve
							unsigned int element = KV_section / (xi_size - 1); // this is the identifier within the bezier element for that NURBS curve. Starts at 0
							unsigned int elem_section = KV_section % (xi_size - 1);    // this is the identifier for the section within that bezier element.  Starts at 0
							
							if (edge_node == 0) {  // update the tri_NURB_elem_section variable but only once per triangle
								vector<unsigned int> tri_NURB_row(5, 0);
								tri_NURB_row[0] = i;
								tri_NURB_row[1] = cur_nurb;
								tri_NURB_row[2] = element;
								tri_NURB_row[3] = elem_section;
								tri_NURB_row[4] = side;
								tri_NURB_elem_section_side.push_back(tri_NURB_row);
							}
							// now that I know that stuff, I need determine what the parametric xi is which corresponds to the physical locations along the edge
							  // first I need to determine the length of the line between the beginning and end of the "half element" 
							vector<double> first_point(3, 0);
							vector<double> second_point(3, 0);
							first_point[0] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section][0];    // this is the x location for the first point of the bounding curve
							first_point[1] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section][1];    // this is the y location for the first point
							first_point[2] = Elem_list[cur_nurb]->elem_Geom[element].weight[elem_section];    // this is the weighting for the first point


							second_point[0] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section + 1][0];    // this is the x location for the second point
							second_point[1] = Elem_list[cur_nurb]->elem_Geom[element].xi_real_loc[elem_section + 1][1];    // this is the y location for the second point
							second_point[2] = Elem_list[cur_nurb]->elem_Geom[element].weight[elem_section + 1];    // this is the weighting for the second point


							double len = sqrt(pow(second_point[0] - first_point[0], 2) + pow(second_point[1] - first_point[1], 2));    // this is the linear length between the beginning and end of the element section


							vector <double> temp(2, 0); // this is a variable I will use to write the value of the current control point location into
							double temp_len = 0;
							int node_index = global_edges[triangles[i]->global_side[side]][edge_node];   // this will let me reference the correct node in the node list

							temp[0] = node_list[node_index][0];
							temp[1] = node_list[node_index][1];
							temp_len = sqrt(pow(temp[0] - first_point[0], 2) + pow(temp[1] - first_point[1], 2));
							temp_len = temp_len / len;
							// now I need to add in the parametric value of the beginning of this section. 
							double begin_xi = Elem_list[cur_nurb]->xi_evals[elem_section];
							double end_xi = Elem_list[cur_nurb]->xi_evals[elem_section + 1];
							temp_len = ((end_xi - begin_xi) * temp_len) + begin_xi;   //  this is the value I need to evaluate the original bezier element at

							temp = eval_Bez_elem(temp_len, element, cur_nurb);

							// now update the weighting field on temp if it is supposed to be first_point or second_point
							double fraction = 1.0 / double(Elem_list[cur_nurb]->p);
							double flag1 = fraction * double(elem_section);
							double flag2 = fraction * double(elem_section + 1);
									
							// determine if the temp point is the same as the first point in the element section
							if (abs(temp_len - flag1) < 0.00001)
								temp.push_back(first_point[2]);
							else if (abs(temp_len - flag2) < 0.00001)
								temp.push_back(second_point[2]);
							else
								temp.push_back(0.0);    // I default to 0 to indicate that I don't know what the weighting is supposed to be

							// now I need to update the position of the node with in node_list
							node_list[node_index] = temp;

							found = true;
						}
						else {   // this means that the Knot vector section is not present within the current NURBS curve
							KV_section -= Elem_list[cur_nurb]->n_el * (xi_size - 1);
							cur_nurb++;
						}
					}
				}

			}
		}

	}
}


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
	// Third, project the coordinates of the nodes
	//
	// Fourth, perfrom curve refinement and most likely degree elevation to reach the final curve
	//
	// Fifth "un-project" the coordinates to get the weights of all of the intermittent points

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

			vector<double> xi_list(2 * degree + 1, 0.0);    // it is this size because there are 2*(degree + 1) - 1 unique nodes
			// first point of the first node is automatically a zero so I don't need to worry about it

			// create vector which will tell the follow loop which nodes to grab, this depends on the side on degree of the triangle
			vector<int> node_side1(degree + 1, 0);
			node_side1[0] = 0;
			node_side1[degree] = 1;
			for (int j = 1; j < degree; j++) {
				node_side1[j] = j + 2;
			}

			vector<int> node_side2(degree, 0);
			node_side2[0] = 1;
			node_side2[degree] = 2;
			for (int j = 1; j < degree; j++) {
				node_side2[j] = j + degree + 1;
			}

			vector<int> node_side3(degree , 0);
			node_side3[0] = 2;
			node_side3[degree] = 0;
			for (int j = 1; j < degree; j++) {
				node_side3[j] = j + (degree << 1);
			}

			// figure out the length of the section, I will do this by linearally connecting the degree + 1 points
			for (unsigned int section = 0; section < same_elem_handle.size(); section++) {
				double section_len = 0.0;
				vector<double> point1(2, 0.0);
				vector<double> point2(2, 0.0);

				for (int j = 0; j < degree; j++) {
					int tri_index = tri_NURB_elem_section_side[same_elem_handle[section]][0];
					//point1[0] = triangles[tri_index]->controlP
				}
			}
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

tri_10_output nurb::tri_10_fast(unsigned int elem, int q)
{

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

