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
#include <numeric>
#include <functional>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


nurb::nurb()
{
	BC.clear();
}


/**********************************************************************************
Function prototype:
void nurb::readSpline(string filename)

Function description:
This function reads in the .spline file which will be specified by its filename
This will read through the file and save the controls, knot vectors, boundary flags,
and boundary condition information into different arrays for later use.

Precondition:
This function will operate correctly if a .spline file is located by the specified
file path.  Furthermore, it must be in the same file format as with Luke Engvall's
MATLAB code.  A full description of this format will be given in the README file

Postcondition:
This sets the value of the point with points to the head of the data structure which
contains the NURBS curve data

**********************************************************************************/


Geom_data* nurb::readSpline(string filename)
{
	// first get the filepath to the current directory
	#ifdef _WIN32
        system("cd > filepath.txt");
    #else
        system("pwd > filepath.txt");
	#endif
	// read that file back in
	ifstream file;
	file.open("filepath.txt");
	getline(file, path_to_file);
	const int end_index = int(path_to_file.size());
	path_to_file.erase(end_index - 5, 5);

	ifstream infile;
	filename += ".spline";
    #ifdef _WIN32
        infile.open(path_to_file + "IO_files\\spline_files\\" + filename);
    #else
        infile.open(path_to_file + "IO_files/spline_files/" + filename);
    #endif

	string line;

	getline(infile, line);  // HEADER line
	getline(infile, line);

	// first important line
	getline(infile, line);
	string temp_num;
	istringstream numStream(line);
	getline(numStream, temp_num, ' ');
	num_curves = stoi(temp_num);            // This value is important to know how many times to loop later

											// Get the first "seed" declaration of the Geom_data structure
	Geom_data *head = new Geom_data;
	Geom_data *current = head;
	current->previous = NULL;
	current->next = NULL;

	getline(infile, line);  //burn a line


							///////////////   Loop through the line for each curve section   //////////////////

	int cur_num;     // this will be the variable of type int that I will reuse when reading in
	for (int i = 0; i < num_curves; i++) {
		/**************************************************************************************************/
		//  First need the data for the number of control points, knot vector length, and polynomial degree
		/**************************************************************************************************/

		getline(infile, line);
		int k = 0;
		int l;
		int j = 0;
		while (j < 3) {
			// loop through all of the white space
			while (line[k] == ' ' || line[k] == '\t') {
				k++;
			}


			// Loop through all the numbers
			l = 0;
			while (line[k + l] != '\t' && line[k + l] != ' ') {
				l++;
			}
			string temp_num = line.substr(k, l);


			cur_num = stoi(temp_num);
			switch (j) {
			case 0:
				current->NCP = cur_num;
				break;
			case 1:
				current->LKV = cur_num;
				break;
			case 2:
				current->p_deg = cur_num;
				break;
			}
			j++;
			k = k + l;
		}

		/*****************************************************/
		// Now go through all of the control point locations
		/*****************************************************/

		getline(infile, line); // burn the 'CONTROL POINTS' line
		rearrange(current, 1, 0);

		for (int j = 0; j < current->NCP; j++) {    // Loop through the rows of the control points array
			getline(infile, line);
			int l = 0;
			int m;
			for (int k = 0; k < 3; k++) {           // Loop through the 3 columns in the spline file

													// loop through all of the white space
				while (line[l] == ' ' || line[l] == '\t') {
					l++;
				}


				// Loop through all the numbers
				m = 0;
				while (line[l + m] != '\t' && line[l + m] != ' ') {
					m++;
					if (l + m == line.length()) {
						break;
					}
				}
				string temp_num = line.substr(l, m);


				double cur_num = atof(temp_num.c_str());
				switch (k) {
				case 0:
					current->controlP[j][k] = cur_num;
					break;

				case 1:
					current->controlP[j][k] = cur_num;
					break;

				case 2:
					current->weight[j] = cur_num;
					break;
				}  // end of the switch statement
				l = l + m;
			}  // end of the k for loop
		}  // end of the j for loop

		   /*****************************************************/
		   //        Now go through the knot vector
		   /*****************************************************/

		getline(infile, line);     // burn the	KNOT VECTOR line
		getline(infile, line);
		l = 0;
		int m;
		for (int j = 0; j < current->LKV; j++) {

			// loop through all of the white space
			while (line[l] == ' ' || line[l] == '\t') {
				l++;
			}


			// Loop through all the numbers
			m = 0;
			while (line[l + m] != '\t' && line[l + m] != ' ') {
				m++;
				if (l + m == line.length()) {
					break;
				}
			}
			string temp_num = line.substr(l, m);

			double cur_num = atof(temp_num.c_str());
			current->KV[j] = cur_num;
			l = l + m;

		}

		// now normalize the knot vector to unit size
		double max_knot = current->KV[current->LKV - 1];
		for (int j = 0; j < current->LKV; j++) {
			current->KV[j] = current->KV[j] / max_knot;
		}

		/*****************************************************/
		//        Now go through the Boundary flags
		/*****************************************************/

		getline(infile, line);    // burn the BOUNDARY FLAGS line
								  // There is no indicator for how long this section will be so I will have to be careful when looping through it
		getline(infile, line);

		while (hasnumber(line)) {   // there are guarenteed to only be two columns so I will just do the loop unrolling myself
			l = 0;
			int m;

			// First number/////////////////////////////////

			// loop through all of the white space
			while (line[l] == ' ' || line[l] == '\t') {
				l++;
			}


			// Loop through all the numbers
			m = 0;
			while (line[l + m] != '\t' && line[l + m] != ' ') {
				m++;
				if (l + m == line.length()) {
					break;
				}
			}
			l = l + m;

			// Second number////////////////////////////////

			// loop through all of the white space
			while (line[l] == ' ' || line[l] == '\t') {
				l++;
			}


			// Loop through all the numbers
			m = 0;
			while (line[l + m] != '\t' && line[l + m] != ' ') {
				m++;
				if (l + m == line.length()) {
					break;
				}
			}
			temp_num = line.substr(l, m);
			int cur_flag = stoi(temp_num);

			// Add this vector into the already created vector
			current->B_flag.push_back(cur_flag);

			getline(infile, line);  // get the next line
		}   // end of the Boundary flags while loop

			/*****************************************************/
			//        Initialize the next Geom_data
			/*****************************************************/
		if (i < num_curves - 1) {
			Geom_data *next_one = new Geom_data;
			next_one->previous = current;
			next_one->next = NULL;
			current->next = next_one;
			current = next_one;   // switch to the next curve
		}


	}
	// Read through the Boundary Conditions section at the end of the file
	getline(infile, line);   // burn the NBOUNDS LINE
	getline(infile, line);

	// read all of the numbers
	int m = 0;
	while (line[m] != '\t' && line[m] != ' ') {
		m++;
		if (m == line.length()) {
			break;
		}
	}
	temp_num = line.substr(0, m);
	int num_bounds = stoi(temp_num);
	rearrange(current, 2, num_bounds);

	// now read through all of the boundary conditions
	for (int i = 0; i < num_bounds; i++) {
		getline(infile, line);
		int l = 0;
		int m;

		// First number/////////////////////////////////

		// loop through all of the white space
		while (line[l] == ' ' || line[l] == '\t') {
			l++;
		}


		// Loop through all the numbers
		m = 0;
		while (line[l + m] != '\t' && line[l + m] != ' ') {
			m++;
			if (l + m == line.length()) {
				break;
			}
		}
		string temp_num = line.substr(l, m);
		int b_num = stoi(temp_num);
		l = l + m;

		// Second number////////////////////////////////

		// loop through all of the white space
		while (line[l] == ' ' || line[l] == '\t') {
			l++;
		}


		// Loop through all the numbers
		m = 0;
		while (line[l + m] != '\t' && line[l + m] != ' ') {
			m++;
			if (l + m == line.length()) {
				break;
			}
		}
		temp_num = line.substr(l, m);
		int b_type = stoi(temp_num);

		BC[i][0] = b_num;
		BC[i][1] = b_type;

	}


	return head;

}

/**********************************************************************************
Function prototype:
void nurb::rearrange(Geom_data *var, string array)

Function description:
This function reshapes the vector array specified by array within the Geom_data
structure var to be the correct width and length.  This code is designed for 2D
problems and therefore only needs 2D arrays

The options for array are
1:     Array of the control points
2:     Array of the boundary conditions


Precondition:
The specified vector array within the Geom_data variable has not been reshaped yet

Postcondition:
The vector arrays within the specified Geom_data variable are reshaped to be the
correct size

**********************************************************************************/

void nurb::rearrange(Geom_data *var, int array, int B_rows)
{
	switch (array) {
	case 1:    // this is the control points
		var->controlP.resize(var->NCP);        // this makes as many rows as the number of control points
		var->weight.resize(var->NCP);          // this makes as many rows as the number of control points but for the weights
		var->KV.resize(var->LKV);              // this makes the knot vector the correct length
		for (int i = 0; i < var->NCP; i++) {   // there are guarenteed to be only 2 columns
			var->controlP[i].resize(2);
		}
	case 2:
		BC.resize(B_rows);        // this makes as many rows as the number of boundary flags
		for (int i = 0; i < B_rows; i++) {    // there are guarenteed to be only 2 columns
			BC[i].resize(2);
		}

	}
}


/**********************************************************************************
Function prototype:
void nurb::hasnumber(string line)

Function description:
This function determines if there is a number present in a given string


Precondition:
The string must be passed in

Postcondition:
This function returns a boolean where true means it contains a number and false
means it does not contain a number

**********************************************************************************/

bool nurb::hasnumber(string line) {
	bool num = false;
	int line_len = int(line.length());
	for (int i = 0; i < line_len; i++) {
		if (isdigit(line[i])) {
			num = true;
		}
	}
	return num;
}


/**********************************************************************************
Function prototype:
void nurb::NURBS2poly(Geom_data *var)

Function description:
This function approximated any number of closed NURBS curves by and equal number of polygons


Precondition:
The head geometric data structure must be passed in

Postcondition:
This function alters the Geom_data so that they contain an array of the nodes of the
polygons and another array describing the connectivity of the polygons

**********************************************************************************/


void nurb::NURBS2poly(Geom_data * var)
{
	while (var != NULL) {
		Bezier_handle *head = extraction1D(var);
		Elem_list.push_back(head);
		vector <double> KV_old = var->KV;
		// now I will update the KV within var to be the C(0) version
		refine_Xi(var);

		// now I need to evaluate the bezier element at p + 1 equispaced points
		evalBez(var, head, -1, KV_old);

		// now determine the arc length of the curve
		for (int i = 0; i < head->n_el; i++) {
			curve_len(head, i);
		}
		// now loop through, dividing up the elements further so that they meet the threshold specifications.
		double thresh = 1.01;      // one percent threshold for the difference between the two curve lengths
		int iter = 1;
		int max_it = 10;
		while (iter <= max_it) {
			bool need_split = false;
			for (int i = 0; i < head->n_el; i++) {
				if (head->elem_Geom[i].curve_len / head->elem_Geom[i].edge_len > thresh) {
					need_split = true;

					// now need to go through curve refinement on that one local element
					// I can "cheat" on the refinement since when you subdivide an element you get two of the same element, only the control points and weights change
					vector<double> xi_to_add(head->p, 0.5);
					curve_refine(head, i, xi_to_add, true);
					i++;    // so I don't subdivide the same element again

				}
				if (head->elem_Geom[i].controlP[0].size() != 3) {   // if no curve refinement was needed, we still need to make it so that the control point have size 3
					for (int m = 0; m <= head->p; m++) {
						head->elem_Geom[i].controlP[m].push_back(head->elem_Geom[i].weight[m]);
					}
				}
			}

			// determine if any splits occured
			if (need_split == true)
				iter++;
			else   // if there weren't any, there is no point to continue looping
				break;

		}

		// now go through and perform a element refine which is based on the geometry of the curve.  This will enhance the resolution where there
		// is high changes in slope or changes in curvature.
		cusp_detection(head);
		shape_refine(head);

		var = var->next;
	}

}


/**********************************************************************************
Function prototype:
void nurb::evalBez(Geom_data *var)

Function description:
This function evaluates a Bezier Element at 10 equispaced locations along the span
0 to 1.

Precondition:
A pointer to the structure Geom_data and Bezier element must be created and passed in

Postcondition:
This function updates various fields of the Bezier_elem struct

**********************************************************************************/

void nurb::evalBez(Geom_data * var, Bezier_handle *Bez, int elem, vector <double> KV_old)
{
	// First I need to determine the knot vector for the element
	int start, stop;

	if (elem == -1) {    // this means the elements need to be created
						 // now I need to evaluate all of the basis functions at p+1 points
		eval_bern(Bez->p, Bez->p + 1, Bez);                                            // need to change this to change where the knot vectors are evaluated

		get_P_and_W(Bez, var, KV_old);

		// now I need to determine the IEN array for the bezier elements in the NURBS curves
		vector<vector<int>>IEN_cur = IEN(var->NCP, var->p_deg, var->KV, Bez->n_el);
		Master_IEN.push_back(IEN_cur);

		// now evaluate each physical location which corresponds to each xi
		stop = Bez->n_el;
		start = 0;
	}
	else {
		stop = elem + 1;
		start = elem;
	}


	for (int i = start; i < stop; i++) {      // loop through all of the elements
		Bez->elem_Geom[i].edge_len = 0.0;
		for (unsigned int j = 0; j < Bez->xi_evals.size(); j++) {      // loop through all of the xi quadrature points
			vector <double> num(2, 0);
			double denom = 0;
			for (int k = 0; k < Bez->p + 1; k++) {            // Loop through all of the control points in that element.  There are p + 1 of them
				num[0] += Bez->N[j][k] * Bez->elem_Geom[i].controlP[k][0] * Bez->elem_Geom[i].weight[k];    // First item is the evaluation of the bernstein polynomial, second is the control point location, third is the weight
				num[1] += Bez->N[j][k] * Bez->elem_Geom[i].controlP[k][1] * Bez->elem_Geom[i].weight[k];
				denom += Bez->N[j][k] * Bez->elem_Geom[i].weight[k];
			}
			num[0] = num[0] / denom;
			num[1] = num[1] / denom;

			vector<double> cur_point;
			cur_point.push_back(num[0]);
			cur_point.push_back(num[1]);
			if (elem == -1)
				Bez->elem_Geom[i].xi_real_loc.push_back(cur_point);     // push this point into the struct
			else
				Bez->elem_Geom[i].xi_real_loc[j] = cur_point;

			// now update the edge length of the NURBS curve
			if (j > 0) {
				Bez->elem_Geom[i].edge_len += sqrt(pow(Bez->elem_Geom[i].xi_real_loc[j - 1][0] - cur_point[0], 2) + pow(Bez->elem_Geom[i].xi_real_loc[j - 1][1] - cur_point[1], 2));      // this is just doing sqrt(dx^2 + dy^2)
			}

		}
	}
}


/**********************************************************************************
Function prototype:
void nurb::curve_len(Bezier_handle *Bez)

Function description:
This function computes the curve length between each adjacent xi value as well as the
linear edge length


Precondition:
Bezier_handle *Bez needs to be passed in.  It also needs the index of the element
for which to evaluate the curve length


Postcondition:
This function updates the curve length field within the struct

**********************************************************************************/

void nurb::curve_len(Bezier_handle *Bez, int i)
{
	int p = Bez->p;
	Bez->elem_Geom[i].curve_len = 0.0;
	int count = 0;
	vector<double> last_point;
	for (int count = 0; count <= (10 * p); count++) {     // this loops through all of the equispaced 10*(p + 1) evaluation points
		double xi = (double)count / (10.0 * (double)p);
		// initialize the numerator and denominator for evaluation
		vector <double> num(2, 0);
		double denom = 0;
		for (int k = 1; k <= p + 1; k++) {           // this loops through all of the basis function
			double N = (fast_fact[p] / (fast_fact[k - 1] * fast_fact[p - k + 1]))*pow(xi, (k - 1)) * pow((1 - xi), (p + 1 - k));
			// the above line is the n choose k representation of Bernstein polynomials over the span [0,1]

			num[0] += N * Bez->elem_Geom[i].controlP[k - 1][0] * Bez->elem_Geom[i].weight[k - 1];    // First item is the evaluation of the bernstein polynomial, second is the control point location, third is the weight
			num[1] += N * Bez->elem_Geom[i].controlP[k - 1][1] * Bez->elem_Geom[i].weight[k - 1];
			denom += N * Bez->elem_Geom[i].weight[k - 1];
		}

		num[0] = num[0] / denom;
		num[1] = num[1] / denom;


		if (count == 0) {
			last_point = num;
		}										// now update the edge length of the NURBS curve
		else {
			Bez->elem_Geom[i].curve_len += sqrt(pow(num[0] - last_point[0], 2) + pow(num[1] - last_point[1], 2));      // this is just doing sqrt(dx^2 + dy^2)
			last_point = num;
		}
	}

}

/**********************************************************************************
Function prototype:
Bezier_elem * nurb::extraction1D(Geom_data *var)

Function description:
This function computes the extraction operators and number of Bezier Elements for
a given NURBS curve


Precondition:
Geom_data *var needs to be passed in.  This is the pointer to the main data struct


Postcondition:
This function updates the extraction operators and number of elements field within
the bezier element struct

**********************************************************************************/
Bezier_handle * nurb::extraction1D(Geom_data *var)
{
	// First initialize the bezier element struct
	Bezier_handle *head = new Bezier_handle;
	MatrixXd Operator;
	Operator.resize(var->p_deg + 1, var->p_deg + 1);
	Operator.setIdentity();


	int a = var->p_deg;
	int b = a + 1;
	int n_el = 1;     // intialize the number of elements
	int n = var->NCP;
	int p = var->p_deg;

	while (b < (n + p)) {    // there is no +1 in here because the indexing starts at 0
							 // first initialize the next operator
		MatrixXd next_one;
		next_one.resize(p + 1, p + 1);
		next_one.setIdentity();


		int i = b;

		while (b < (n + p) && (var->KV[b + 1] == var->KV[b])) {
			b++;
		}

		int mult = b - i + 1;    //compute the multiplicity of b
		if (mult < p) {
			double numer = var->KV[b] - var->KV[a];

			// create vector alpha_s
			vector <double> alpha_s(p - mult);
			for (int j = p; j >= mult + 1; j--) {
				alpha_s[j - mult - 1] = numer / (var->KV[a + j] - var->KV[a]);
			}

			int r = p - mult;     // compute the number of knots needed to be added

			for (int j = 1; j <= r; j++) {
				int save = r - j + 1;
				int s = mult + j;

				for (int k = p; k >= s; k--) {
					double alpha = alpha_s[k - s];
					Operator.col(k) = alpha*Operator.col(k) + (1 - alpha)* Operator.col(k - 1);
				}

				if (b < (n + p)) {
					for (int i = 0; i <= j; i++) {
						next_one(save + i - 1, save - 1) = Operator(p - j + i, p);
					}
				}
			}
		}

		// setup sext element
		if (b < (n + p)) {
			a = b;
			b = a + 1;
			n_el++;
		}
		head->Operator.push_back(Operator);
		Operator = next_one;
	}


	head->n_el = n_el;
	head->p = var->p_deg;
	return head;
}


/**********************************************************************************
Function prototype:
void nurb::fact_table()

Function description:
This function creates a table of common factorials, UP TO 15!


Precondition:
This function is called at the beginning of mesh_main.cpp and the data is saved for
the remainder of the program's execution

Postcondition:
This function saves the table within the protected fields of the class
**********************************************************************************/
void nurb::fact_table(int p)
{
	double total = 1;
	fast_fact.push_back(total); // this is 0!
	for (double i = 1; i <= p; i++) {
		total *= i;
		fast_fact.push_back(total);
	}
}


/**********************************************************************************
Function prototype:
void nurb::eval_bern(int p, int n, Bezier_elem *Bez)

Function description:
This function evaluates the bernstein basis functions at a set of equispaced
points.  There are p+1 of these points


Precondition:
This function requires the polynomial degree of the element (p), the number of
basis functions in the element (n), and a point to its parent Bezier_elem struct

Postcondition:
This function updates the bezier_value struct within the parent Bezier_elem
**********************************************************************************/

void nurb::eval_bern(int p, int n, Bezier_handle *Bez)
{
	vector <double> xi_list;

	for (int count = 0; count <= p; count++) {     // this loops through all of the equispaced p + 1 evaluation points  I CHANGED THIS FROM 2 to p
		double xi = (double)count / p;            // this is so it evaluates equispace p + 1 points
		xi_list.push_back(xi);
		vector <double> same_xi;
		for (int i = 1; i <= p + 1; i++) {           // this loops through all of the basis function
			double N = (fast_fact[p] / (fast_fact[i - 1] * fast_fact[p - i + 1]))*pow(xi, (i - 1)) * pow((1 - xi), (p + 1 - i));
			// the above line is the n choose k representation of Bernstein polynomials over the span [0,1]
			same_xi.push_back(N);
		}
		Bez->N.push_back(same_xi);
	}
	Bez->xi_evals = xi_list;
}

/**********************************************************************************
Function prototype:
void nurb::get_P_and_W(Bezier_elem *Bez, Geom_data *var)

Function description:
This function determines the control points and weights which are associated with
each bezier element


Precondition:
This function requires a pointer to the Geom_data struct, which houses the original
geometric data, as well as a pointer to the current Bezier Element, Bez.

Postcondition:
This function updates the bezier_value struct within the parent Bezier_elem and
return a boolean.  If true, it means that the program had to create all of the data
like it would for the first time this function is called.  If it is false, it means
that the function didn't do anything since there is already data.
**********************************************************************************/

void nurb::get_P_and_W(Bezier_handle * Bez, Geom_data *var, vector<double> KV_old)
{
	// First pull off the current IEN array
	int n = int(KV_old.size()) - (Bez->p + 1);
	vector <vector <int>> IEN_cur = IEN(n, Bez->p, KV_old, Bez->n_el);

	// Next construct the projected control points into the matrix data type
	for (int i = 0; i < Bez->n_el; i++) {
		// I need to make a seperate matrix for the X Y and Z coordinate.
		MatrixXd proj;   // first column is x, second column is y, third is z (weight)
		proj.resize(Bez->p + 1, 3);

		// Assemble the matrix
		for (int j = 0; j <= Bez->p; j++) {
			for (int k = 0; k < 3; k++) {
				if (k != 2) {
					proj(j, k) = var->controlP[IEN_cur[i][j]][k];
				}
				else {
					proj(j, k) = var->weight[IEN_cur[i][j]];
				}
			}
		}
		// now need to multiply this vector by the transpose of the extraction operator
		MatrixXd Bez_points;

		Bez_points.resize(Bez->p + 1, 3);

		Bez_points = Bez->Operator[i].transpose()*proj;

		Bez_points.transposeInPlace();     // this is so it is essentially in row major order instead of column major

										   // Now I extract the x,y, and weight data and store it appropriately
		Bezier_elem proj_data;
		double temp;
		for (int j = 0; j < Bez_points.cols(); j++) {
			vector<double> cur_point;
			for (int k = 0; k < 3; k++) {
				temp = *(Bez_points.data() + (k + (j * 3)));    // get the value from that cell in the matrix

				if (k != 2) {   // if it is supposed to be either the x or y coordinate
					cur_point.push_back(temp);
				}
				else
					proj_data.weight.push_back(temp);
			}
			proj_data.controlP.push_back(cur_point);
			proj_data.side = i;
		}
		Bez->elem_Geom.push_back(proj_data);  // push all of this data onto the vector
	}
}

/**********************************************************************************
Function prototype:
void nurb::IEN(int n; int p; vector<double> Xi)

Function description:
This function computes the IEN which describes the local to global connectivity
of the Bezier control points/basis functions.


Precondition:
This function requires the number of basis functions in the C(0) NURBS curve,
the degree for this same NURBS curve, and its C(0) Knot vector

Postcondition:
This function updates the Master_IEN protected data variable
**********************************************************************************/

vector<vector<int>> nurb::IEN(int n, int p, vector<double> Xi, int n_el)
{
	// this will make the transpose of the Matlab version of this code.  The reason this is done is because c++ has row order to it while Matlab has column order
	// the numbers within the IEN array denote the index of the basis function, therefore, it will start at 0 instead of 1
	vector <vector<int>> IEN(n_el);
	int l = p + 1;
	int e = 0;
	while (l <= n) {
		for (int a = 0; a < p + 1; a++) {
			IEN[e].push_back((l + a) - (p + 1));
		}
		l++;
		while (Xi[l] == Xi[l - 1] && (l < n)) {
			l++;
		}
		if (l <= n)
			e++;
	}
	// now combine into the master array with stores the IEN Arrays for each of the NURBS curves
	return IEN;
}

/**********************************************************************************
Function prototype:
void nurb::subdivide_element(Bezier_handle *Bez, int e_index)

Function description:
This function divides the given element in half, creates is corresponding Extraction
Operator from the original NURBS curve and updates all appropriate data field within
the Bezier_handler struct


Precondition:
This function requires a pointer to the Bezier_handler struct as well as an index
corresponding to the element which is to be subdivided

Postcondition:
This function updates the Bezier_handler *Bez
**********************************************************************************/

void nurb::subdivide_element(Bezier_handle * Bez, int e_index)
{
	// first create the knot vector
	vector<double> Xi;
	for (int i = 0; i < Bez->p + 1; i++)
		Xi.push_back(0.0);
	for (int i = 0; i < Bez->p + 1; i++)
		Xi.push_back(1.0);
	// now perform a curve refine to get the control points and weights for the refined version
	for (int added_points = 0; added_points < Bez->p; added_points++) {  // loop through all of the points that need to be added
																		 // first determine the length variable m

		int m = int(Bez->elem_Geom[e_index].controlP.size()) + 1;    // m is the number of control points after I add one the new one
		int k = Bez->p + 1;      // this is the location of the added point
		vector<vector<double>> Q(m, vector<double>(2, 0.0));
		vector<double> W(m, 0);

		// first I need to project the curve
		for (int j = 1; j <= m; j++) {
			if (j == 1) {
				Q[j - 1][0] = Bez->elem_Geom[e_index].controlP[0][0] * Bez->elem_Geom[e_index].weight[j - 1];     // I am multiplying by the weights to project the curve
				Q[j - 1][1] = Bez->elem_Geom[e_index].controlP[0][1] * Bez->elem_Geom[e_index].weight[j - 1];
				W[j - 1] = Bez->elem_Geom[e_index].weight[m - 2];

			}
			else if (j == m) {
				Q[j - 1][0] = Bez->elem_Geom[e_index].controlP[m - 2][0] * Bez->elem_Geom[e_index].weight[m - 2];
				Q[j - 1][1] = Bez->elem_Geom[e_index].controlP[m - 2][1] * Bez->elem_Geom[e_index].weight[m - 2];
				W[j - 1] = Bez->elem_Geom[e_index].weight[m - 2];

			}
			else {
				double alpha;
				if (j >= 1 && j <= k - Bez->p) {
					alpha = 1.0;
				}
				else if (j >= k - Bez->p + 1 && j <= k) {
					alpha = (0.5 - Xi[j - 1]) / (Xi[j + Bez->p - 1] - Xi[j - 1]);
				}
				else
					alpha = 0.0;
				// now save the values of the data points and project them
				double point11 = Bez->elem_Geom[e_index].controlP[j - 1][0] * Bez->elem_Geom[e_index].weight[j - 1];   //11 = current j index, x or y (x = 0, y = 1)
				double point12 = Bez->elem_Geom[e_index].controlP[j - 1][1] * Bez->elem_Geom[e_index].weight[j - 1];
				double point21 = Bez->elem_Geom[e_index].controlP[j - 2][0] * Bez->elem_Geom[e_index].weight[j - 2];
				double point22 = Bez->elem_Geom[e_index].controlP[j - 2][1] * Bez->elem_Geom[e_index].weight[j - 2];

				Q[j - 1][0] = alpha*point11 + (1 - alpha)*point21;
				Q[j - 1][1] = alpha*point12 + (1 - alpha)*point22;
				W[j - 1] = alpha*Bez->elem_Geom[e_index].weight[j - 1] + (1 - alpha)*Bez->elem_Geom[e_index].weight[j - 2];

			}
			// now come back from projected space

			Q[j - 1][0] = Q[j - 1][0] / W[j - 1];
			Q[j - 1][1] = Q[j - 1][1] / W[j - 1];

		}



		Bez->elem_Geom[e_index].controlP = Q;
		Bez->elem_Geom[e_index].weight = W;
		Xi.insert(Xi.begin() + Bez->p + 1, 0.5);

	}

	// now I need to seperate the 1 element into 2
	// the knot vector has C(0) continuity at this point so I don't have to do a full extraction
	// also, therefore the extraction operator between the 2 element version and the one element is
	// the identity matrix when multiplying the old and new Extraction operators together to obtain a
	// NURBS to refined element operator, you would be multiplying the other matrix by the identity matrix
	// therefore, that step is largely ignored.

	Bez->n_el++;


	// so simply make a copy of the previous Operator
	Bez->Operator.insert(Bez->Operator.begin() + e_index, Bez->Operator[e_index]);
	Bezier_elem inserted = Bez->elem_Geom[e_index];      // first just copy the element to be divided
	Bez->elem_Geom.insert(Bez->elem_Geom.begin() + e_index, inserted);    // this makes the inserted element come before the one that was determined to be repeated

																		  // Now Update the master IEN array so it reflects the added control points.
	vector<int> insert(Bez->p + 1, 0);  // this is the entry corresponding to the inserted point
	vector<vector<int>> cur_IEN = Master_IEN.back();
	insert[0] = cur_IEN[e_index][Bez->p];
	for (int i = 1; i < Bez->p + 1; i++) {
		insert[i] = insert[i - 1] + 1;
	}
	cur_IEN.insert(cur_IEN.begin() + e_index + 1, insert);


	for (int i = e_index + 2; i < Bez->n_el; i++) {
		for (int j = 0; j < Bez->p + 1; j++) {
			cur_IEN[i][j] += Bez->p;
		}

	}
	Master_IEN.back() = cur_IEN;

	// now update the weights and control points on the new one
	// for this element, I will just pop the last p elements off the back
	for (int i = 0; i < Bez->p; i++) {
		Bez->elem_Geom[e_index].controlP.pop_back();
		Bez->elem_Geom[e_index].weight.pop_back();
	}

	// now update the weights and control points for the old one.
	// this will contain control points that occur later.
	// this will remove the elements from the front
	for (int i = 0; i < Bez->p; i++) {
		Bez->elem_Geom[e_index + 1].controlP.erase(Bez->elem_Geom[e_index + 1].controlP.begin());
		Bez->elem_Geom[e_index + 1].weight.erase(Bez->elem_Geom[e_index + 1].weight.begin());
	}


	// now I need to find the curve and edge lengths of these new elements
	Geom_data *garbage = new Geom_data;    // just need one of these so it will be happy
	vector<double> trash(1, 0);
	evalBez(garbage, Bez, e_index, trash);
	evalBez(garbage, Bez, e_index + 1, trash);
	curve_len(Bez, e_index);
	curve_len(Bez, e_index + 1);
	delete garbage;
}

/**********************************************************************************
Function prototype:
void nurb::refine_Xi(Geom_data *var)

Function description:
This function alters the knot vector within var so that it exihibits C(0)
continuity


Precondition:
This function only requires a valid pointer to Geom_data struct in order to perform

Postcondition:
This function updates var->KV
**********************************************************************************/

void nurb::refine_Xi(Geom_data * var)
{
	vector<double> Xi = var->KV;
	unsigned int i = 0;
	while (Xi[i] == Xi[0]) {
		i++;
	}
	while (true) {
		int start = i;
		// loop through to find the end of this same reoccuring number
		while (Xi[i] == Xi[start] && i < Xi.size() - 1) {
			i++;
		}
		unsigned int amount_cur = i - start;
		unsigned int need = var->p_deg - amount_cur;
		if (i == Xi.size() - 1) {
			var->KV = Xi;
			var->LKV = int(Xi.size());
			var->NCP = var->LKV - (var->p_deg + 1);
			return;
		}
		else {
			vector <double> insert_num(need, Xi[start]);
			Xi.insert(Xi.begin() + start, insert_num.begin(), insert_num.end());
			i += need;
		}
	}
}

/**********************************************************************************
Function prototype:
void nurb::shape_refine(Bezier_handle * Bez)

Function description:
This function evaluates the curvature at xi values of 0.2, 0.4, 0.6, and 0.8 for every
Bezier element along the NURBS curve and compiles it into a list of increasing curvature.
Then it goes through that list and if any of these locations has a curvature which is
10 * (the median) + 1, curve refinement is implimented to divide the Bezier element.

**********************************************************************************/

void nurb::shape_refine(Bezier_handle * Bez)
{
	vector<pair<double, int>> kurv_holder;
	int count = 0;
	for (int i = 0; i < Bez->n_el; i++) {    // loop through all of the elements in the curve
											 // now find the curvature at 4 equispaced points within the element
		for (int j = 0; j < 4; j++) {
			double kurv = get_kurv(Bez, i, double(j + 1) / 5);
			kurv_holder.push_back(make_pair(kurv, count));
			count++;
		}
	}
	// now find the median value
	sort(kurv_holder.begin(), kurv_holder.end());
	int index = (int(kurv_holder.size()) / 2) - 1;

	// now make a list of all of the indexes that need splitting
	vector<int> split_index;
	for (int i = index + 1; i < int(kurv_holder.size()); i++) {
		if (kurv_holder[i].first >(10 * kurv_holder[index].first) + 1)
			split_index.push_back(kurv_holder[i].second);
	}
	sort(split_index.begin(), split_index.end());
	int cur_elem;
	vector<double> xi_to_add;
	vector<double> xi_choices(4, 0);
	xi_choices[0] = 0.2;
	xi_choices[1] = 0.4;
	xi_choices[2] = 0.6;
	xi_choices[3] = 0.8;


	for (int i = 0; i < int(split_index.size()); i++) {
		if (xi_to_add.size() == 0) { // there is nothing in here so I can just add the next one
			cur_elem = split_index[i] / 4;
			double val = xi_choices[split_index[i] % 4];
			for (int j = 0; j < Bez->p; j++)
				xi_to_add.push_back(val);
		}
		else {
			if ((split_index[i] / 4) != cur_elem) {  // the xi_to_add is full so I need to split up the Bezier element
				curve_refine(Bez, cur_elem, xi_to_add, true);
				int num_add = int(xi_to_add.size());
				// now I need to update the rest of split list
				for (int j = i; j < int(split_index.size()); j++) {
					split_index[j] += 4 * (num_add / Bez->p);

				}
				xi_to_add.clear();
				i--;
			}
			else {  // this means that we are still in the current element so I can just add the next point in
				double val = xi_choices[split_index[i] % 4];
				for (int j = 0; j < Bez->p; j++)
					xi_to_add.push_back(val);
			}
		}

	}
	if (xi_to_add.size())  // handle the last case if there are ones to do
		curve_refine(Bez, cur_elem, xi_to_add, true);   // take care of the last one


}


/**********************************************************************************
Function prototype:
void nurb::get_kurv(Bezier_handle * Bez, int elem, double xi)

Function description:
This function calculates the curvature of the specified NURBS curve, Bezier element,
and xi value and returns it to shape_refine.  It follows the process and variables from 
M. S. Floater in this paper http://www.mn.uio.no/math/english/people/aca/michaelf/papers/bez.pdf
The section of importantance is Corollary 6 on page 14

**********************************************************************************/
double nurb::get_kurv(Bezier_handle * Bez, int elem, double xi)
{
	if (Bez->p < 2)
		return 0;
	// start with figuring out the n-2 ones
	int k = Bez->p - 2;
	vector<double> p0n2(2, 0);
	vector<double> p1n2(2, 0);
	vector<double> p2n2(2, 0);
	double w0n2 = 0.0;
	double w1n2 = 0.0;
	double w2n2 = 0.0;

	for (int j = 0; j <= k; j++) {
		double B = n_choose_k(k, j) * pow(xi, j) * pow(1 - xi, k - j);
		p0n2[0] = p0n2[0] + B * Bez->elem_Geom[elem].controlP[j][0] * Bez->elem_Geom[elem].controlP[j][2];
		p0n2[1] = p0n2[1] + B * Bez->elem_Geom[elem].controlP[j][1] * Bez->elem_Geom[elem].controlP[j][2];
		p1n2[0] = p1n2[0] + B * Bez->elem_Geom[elem].controlP[j + 1][0] * Bez->elem_Geom[elem].controlP[j + 1][2];
		p1n2[1] = p1n2[1] + B * Bez->elem_Geom[elem].controlP[j + 1][1] * Bez->elem_Geom[elem].controlP[j + 1][2];
		p2n2[0] = p2n2[0] + B * Bez->elem_Geom[elem].controlP[j + 2][0] * Bez->elem_Geom[elem].controlP[j + 2][2];
		p2n2[1] = p2n2[1] + B * Bez->elem_Geom[elem].controlP[j + 2][1] * Bez->elem_Geom[elem].controlP[j + 2][2];

		w0n2 += B * Bez->elem_Geom[elem].controlP[j][2];
		w1n2 += B * Bez->elem_Geom[elem].controlP[j + 1][2];
		w2n2 += B * Bez->elem_Geom[elem].controlP[j + 2][2];
	}
	p0n2[0] = p0n2[0] / w0n2;
	p0n2[1] = p0n2[1] / w0n2;

	p1n2[0] = p1n2[0] / w1n2;
	p1n2[1] = p1n2[1] / w1n2;

	p2n2[0] = p2n2[0] / w2n2;
	p2n2[1] = p2n2[1] / w2n2;

	// now do all of the n-1 ones
	k++;
	vector<double> p0n1(2, 0);
	vector<double> p1n1(2, 0);
	double w0n1 = 0.0;
	double w1n1 = 0.0;

	for (int j = 0; j <= k; j++) {
		double B = n_choose_k(k, j) * pow(xi, j) * pow(1 - xi, k - j);
		p0n1[0] += B * Bez->elem_Geom[elem].controlP[j][0] * Bez->elem_Geom[elem].controlP[j][2];
		p0n1[1] += B * Bez->elem_Geom[elem].controlP[j][1] * Bez->elem_Geom[elem].controlP[j][2];
		p1n1[0] += B * Bez->elem_Geom[elem].controlP[j + 1][0] * Bez->elem_Geom[elem].controlP[j + 1][2];
		p1n1[1] += B * Bez->elem_Geom[elem].controlP[j + 1][1] * Bez->elem_Geom[elem].controlP[j + 1][2];

		w0n1 += B * Bez->elem_Geom[elem].controlP[j][2];
		w1n1 += B * Bez->elem_Geom[elem].controlP[j + 1][2];
	}
	p0n1[0] = p0n1[0] / w0n1;
	p0n1[1] = p0n1[1] / w0n1;
	p1n1[0] = p1n1[0] / w1n1;
	p1n1[1] = p1n1[1] / w1n1;


	// now do the n ones
	k++;
	double w0n0 = 0.0;
	for (int j = 0; j <= k; j++) {
		double B = n_choose_k(k, j) * pow(xi, j) * pow(1 - xi, k - j);
		w0n0 += B * Bez->elem_Geom[elem].controlP[j][2];
	}

	// now prepare some variables to make the curvature formula easier
	vector<double> left_cross(2, 0);
	vector<double> right_cross(2, 0);
	vector<double> bottom(2, 0);

	left_cross[0] = p1n2[0] - p0n2[0];
	left_cross[1] = p1n2[1] - p0n2[1];

	right_cross[0] = p2n2[0] - p1n2[0];
	right_cross[1] = p2n2[1] - p1n2[1];

	bottom[0] = p1n1[0] - p0n1[0];
	bottom[1] = p1n1[1] - p0n1[1];

	double R1 = (w0n2 * w1n2 * w2n2 * pow(w0n0, 3)) / (pow(w0n1, 3) * pow(w1n1, 3));
	double cross = abs((left_cross[0] * right_cross[1]) - (left_cross[1] * right_cross[0]));
	double denom = pow(sqrt(pow(bottom[0], 2) + pow(bottom[1], 2)), 3);
	double kurv = (double(k - 1) / double(k)) * R1 * cross / denom;
	return kurv;

}

/**********************************************************************************
Function prototype:
void nurb::cusp_detection(Bezier_handle * Bez)

Function description:
This function determines first derivative vectors at xi = 0.001 and xi = 0.999 in adjacent
elements so that these two location are close to each other.  Then upon employing the
dot product formula, a turning angle is found between these two vector.  Then both of these
two element are refine according to the following turning angle scheme (numbers listed
are in degrees).

100 <  angle < 130    =    1 equispaced refinement point each element
130 <= angle < 145    =    2 equispaced refinement point each element
145 <= angle < 155    =    3 equispaced refinement point each element
155 <= angle < 165    =    4 equispaced refinement point each element
165 <= angle < 175    =    5 equispaced refinement point each element
175 <= angle < 180    =    6 equispaced refinement point each element

**********************************************************************************/
void nurb::cusp_detection(Bezier_handle * Bez)
{
	// I am going to do this by getting first derivative vectors at xi = 0.01 and 0.99 for each element
	// and then determining the turning angle between each inter-element pair.
	// Then I will have the number of inserted knots depending on where it is between 100 and 180 degrees
	vector<vector<double>> derivs;
	for (int i = 0; i < Bez->n_el; i++) {
		derivs.push_back(get_deriv(Bez, i, 0.001));
		derivs.push_back(get_deriv(Bez, i, 0.999));
	}
	// now I need to move the 0th entry to be the last entry so that it is next to its pair
	derivs.push_back(derivs[0]);
	// now I need to erase the 0th entry from the front of the list
	derivs.erase(derivs.begin());

	// now I need to compute the turning angle
	vector<double> turning;
	for (unsigned int i = 0; i < derivs.size(); i++) {
		double dot = (derivs[i][0] * derivs[i + 1][0]) + (derivs[i][1] * derivs[i + 1][1]);
		// I don't need the denominator because I already know that they are unit sized
		turning.push_back(acos(dot) * 180 / 3.141592654);
		i++;
	}

	// now I need to go back through the turning vector and refine the curve depending on the size of its turning angle
	for (unsigned int i = 0; i < turning.size(); i++) {
		if (turning[i] > 100) {  // then I am calling it a cusp
			vector<double> xi_to_add;

			int second_index = i + 1;
			if (i == turning.size() - 1)   // reconcile that fact that I moved the first entry in turning to be the last one
				second_index = 0;

			if (turning[i] < 130) {  // it is between 100 and 130
				for (int j = 0; j < Bez->p; j++)
					xi_to_add.push_back(0.5);
				curve_refine(Bez, i, xi_to_add, true);
				if (second_index)
					second_index++;
				curve_refine(Bez, second_index, xi_to_add, true);
			}
			else if (turning[i] < 145) {  // it is between 130 and 145
				for (int j = 0; j <= Bez->p; j++) {
					xi_to_add.push_back(1.0 / 3.0);
					xi_to_add.push_back(2.0 / 3.0);
				}
				sort(xi_to_add.begin(), xi_to_add.end());
				curve_refine(Bez, i, xi_to_add, true);
				if (second_index)
					second_index += 2;
				curve_refine(Bez, second_index, xi_to_add, true);
			}
			else if (turning[i] < 155) {  // it is between 145 and 155
				for (int j = 0; j < Bez->p; j++) {
					xi_to_add.push_back(1.0 / 4.0);
					xi_to_add.push_back(2.0 / 4.0);
					xi_to_add.push_back(3.0 / 4.0);
				}
				sort(xi_to_add.begin(), xi_to_add.end());
				curve_refine(Bez, i, xi_to_add, true);
				if (second_index)
					second_index += 3;
				curve_refine(Bez, second_index, xi_to_add, true);
			}
			else if (turning[i] < 165) {  // it is between 155 and 165
				for (int j = 0; j < Bez->p; j++) {
					xi_to_add.push_back(1.0 / 5.0);
					xi_to_add.push_back(2.0 / 5.0);
					xi_to_add.push_back(3.0 / 5.0);
					xi_to_add.push_back(4.0 / 5.0);
				}
				sort(xi_to_add.begin(), xi_to_add.end());
				curve_refine(Bez, i, xi_to_add, true);
				if (second_index)
					second_index += 4;
				curve_refine(Bez, second_index, xi_to_add, true);
			}
			else if (turning[i] < 175) {  // it is between 165 and 175
				for (int j = 0; j < Bez->p; j++) {
					xi_to_add.push_back(1.0 / 6.0);
					xi_to_add.push_back(2.0 / 6.0);
					xi_to_add.push_back(3.0 / 6.0);
					xi_to_add.push_back(4.0 / 6.0);
					xi_to_add.push_back(5.0 / 6.0);
				}
				sort(xi_to_add.begin(), xi_to_add.end());
				curve_refine(Bez, i, xi_to_add, true);
				if (second_index)
					second_index += 5;
				curve_refine(Bez, second_index + 5, xi_to_add, true);
			}
			else if (turning[i] < 180) {   // this could essentially be an else statement as the angle will never be more than 180
				for (int j = 0; j < Bez->p; j++) {
					xi_to_add.push_back(1.0 / 7.0);
					xi_to_add.push_back(2.0 / 7.0);
					xi_to_add.push_back(3.0 / 7.0);
					xi_to_add.push_back(4.0 / 7.0);
					xi_to_add.push_back(5.0 / 7.0);
					xi_to_add.push_back(6.0 / 7.0);
				}
				sort(xi_to_add.begin(), xi_to_add.end());
				curve_refine(Bez, i, xi_to_add, true);
				if (second_index)
					second_index += 6;
				curve_refine(Bez, second_index, xi_to_add, true);

			}
		}
	}
}

/**********************************************************************************
Function prototype:
void nurb::get_deriv(Bezier_handle * Bez, int elem, double xi)

Function description:
This function is the helper function to cusp_detection which determines the 2
coordinate first derivative vectors which will then be used in the dot product
formula.  This also follows the structure of M. S. Floater formula in
http://www.mn.uio.no/math/english/people/aca/michaelf/papers/bez.pdf

**********************************************************************************/
vector<double> nurb::get_deriv(Bezier_handle * Bez, int elem, double xi)
{
	vector<double> result(2, 0);

	// start with figuring out the n-1 ones
	int k = Bez->p - 1;
	vector<double> p0n1(2, 0);
	vector<double> p1n1(2, 0);
	double w0n1 = 0.0;
	double w1n1 = 0.0;

	for (int j = 0; j <= k; j++) {
		double B = n_choose_k(k, j) * pow(xi, j) * pow(1 - xi, k - j);
		p0n1[0] += B * Bez->elem_Geom[elem].controlP[j][0] * Bez->elem_Geom[elem].controlP[j][2];
		p0n1[1] += B * Bez->elem_Geom[elem].controlP[j][1] * Bez->elem_Geom[elem].controlP[j][2];
		p1n1[0] += B * Bez->elem_Geom[elem].controlP[j + 1][0] * Bez->elem_Geom[elem].controlP[j + 1][2];
		p1n1[1] += B * Bez->elem_Geom[elem].controlP[j + 1][1] * Bez->elem_Geom[elem].controlP[j + 1][2];

		w0n1 += B * Bez->elem_Geom[elem].controlP[j][2];
		w1n1 += B * Bez->elem_Geom[elem].controlP[j + 1][2];
	}
	p0n1[0] = p0n1[0] / w0n1;
	p0n1[1] = p0n1[1] / w0n1;
	p1n1[0] = p1n1[0] / w1n1;
	p1n1[1] = p1n1[1] / w1n1;


	// now do the n ones
	k++;
	double w0n0 = 0.0;
	for (int j = 0; j <= k; j++) {
		double B = n_choose_k(k, j) * pow(xi, j) * pow(1 - xi, k - j);
		w0n0 += B * Bez->elem_Geom[elem].controlP[j][2];
	}

	result[0] = k * (w0n1 * w1n1 / pow(w0n0, 2)) * (p1n1[0] - p0n1[0]);
	result[1] = k * (w0n1 * w1n1 / pow(w0n0, 2)) * (p1n1[1] - p0n1[1]);

	double abs_result = sqrt(pow(result[0], 2) + pow(result[1], 2));

	result[0] = result[0] / abs_result;
	result[1] = result[1] / abs_result;
	return result;
}
