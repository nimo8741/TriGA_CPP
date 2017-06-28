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

void nurb::smoothMesh()
{
	// first I need to find all of the boundary triangles and adjust all of the boundary control points so the line exactly on the geometry
	adjust_boundary();
	// now I need to solve the linear elasticity
	LE2D();
}

void nurb::adjust_boundary()
{
	for (unsigned int i = 0; i < triangles.size(); i++) {  //  loop through all of the triangles
		for (int side = 0; side < 3; side++) {
			if (global_edges[triangles[i]->global_side[side]][4] != phy_groups) {  // if it passes this, it is a boundary edge

																				   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
																				   // add the nodes on this edge to the bNodes array if they are not repeats
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

void nurb::LE2D()
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
	evaluate_tri_basis();

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
			tri_10_output current = tri_10_fast(q);

		}

	}
}

void nurb::evaluate_tri_basis()
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
	MatrixXi tuples(10, 3);
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

tri_10_output nurb::tri_10_fast(int q)
{
	unsigned int nen = tri_N[0].size();
	tri_10_output output;   // create instance of struct I will eventually return
	output.R.resize(nen);

	return output;
}

