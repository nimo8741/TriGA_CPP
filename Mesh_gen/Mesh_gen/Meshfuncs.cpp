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
	ifstream infile;
	infile.open(filename);
	string line;

	getline(infile, line);  // HEADER line
	getline(infile, line);

	// first important line
	getline(infile, line);
	string temp_num;
	istringstream numStream(line);
	getline(numStream, temp_num,' ');
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
				k = k++;
			}


			// Loop through all the numbers
			l = 0;
			while (line[k + l] != '\t' && line[k+l] != ' ') {
				l = l++;
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
					l = l++;
				}


				// Loop through all the numbers
				m = 0;
				while (line[l + m] != '\t' && line[l + m] != ' ') {
					m = m++;
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
				l = l++;
			}


			// Loop through all the numbers
			m = 0;
			while (line[l + m] != '\t' && line[l + m] != ' ') {
				m = m++;
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
				l = l++;
			}


			// Loop through all the numbers
			m = 0;
			while (line[l + m] != '\t' && line[l + m] != ' ') {
				m = m++;
				if (l + m == line.length()) {
					break;
				}
			}
			string temp_num = line.substr(l, m);
			int cur_side = stoi(temp_num);
			l = l + m;

			// Second number////////////////////////////////

			// loop through all of the white space
			while (line[l] == ' ' || line[l] == '\t') {
				l = l++;
			}


			// Loop through all the numbers
			m = 0;
			while (line[l + m] != '\t' && line[l + m] != ' ') {
				m = m++;
				if (l + m == line.length()) {
					break;
				}
			}
			temp_num = line.substr(l, m);
			int cur_flag = stoi(temp_num);

			// Combine both of these numbers into a small vector
			vector<int> new_row;
			new_row.push_back(cur_side);
			new_row.push_back(cur_flag);

			// Add this vector into the already created vector
			current->B_flag.push_back(new_row);
			
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
		m = m++;
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
			l = l++;
		}


		// Loop through all the numbers
		m = 0;
		while (line[l + m] != '\t' && line[l + m] != ' ') {
			m = m++;
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
			l = l++;
		}


		// Loop through all the numbers
		m = 0;
		while (line[l + m] != '\t' && line[l + m] != ' ') {
			m = m++;
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
	int line_len = line.length();
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
		// First create the xi vector
		int knot_len = var->KV.size();

		var->xi.push_back(var->KV[0]);
		for (int i = 0; i < knot_len; i++) {
			if (var->KV[i] != var->xi.back()) {
				var->xi.push_back(var->KV[i]);
			}
		}
		// now add in the intermediate values
		knot_len = var->xi.size();
		for (int i = 0; i < knot_len - 1; i++) {
			var->xi.insert(var->xi.begin()+(2*i +1), 0.5*(var->xi[2*i] + var->xi[2*i+1]));
		}

		Bezier_elem *head = extraction1D(var);
		Elem_list.push_back(head);
		evalNURBS(var);
		//curve_len(var);
		var = var->next;
	}


}



/**********************************************************************************
Function prototype:
void nurb::evalNURBS(Geom_data *var)

Function description:
This function evaluates a Bezier Element at 10 equispaced locations along the span
0 to 1.

Precondition:
A pointer to the structure Geom_data had to be created and passed in

Postcondition:
This function updates the node field within the Geom_data with the physical 
coordinate of each xi evaluation point

**********************************************************************************/

void nurb::evalNURBS(Geom_data * var)
{
	make_family(var->p_deg, var->NCP - 1, var->KV, var->xi);  // the reason there is the minus 1 is since the indexing begins at 0 instead of 1

	// now traverse the structure to evaluate each physical location which corresponds to each xi
	vector<deBoor_fam *> curList = deBoor_tree.back();
	int len = curList.size();
	for (int i = 0; i < len; i++) {
		vector<double> num(2, 0.0);   // there are two values here since there is both the x and y location
		double denom = 0.0;           // this is just a single scalar
		deBoor_fam *current = curList[i];
		int count = 0;
		while (current != NULL){
			num[0] = num[0] + current->N * var->weight[count] * var->controlP[count][0];
			num[1] = num[1] + current->N * var->weight[count] * var->controlP[count][1];
			denom = denom + current->N * var->weight[count];
			count++;
			current = current->brother;
		}
		num[0] = num[0] / denom;   // I will reuse this variable to make it slightly simpler
		num[1] = num[1] / denom;
		var->node.push_back(num);
	}
}

/**********************************************************************************
Function prototype:
void nurb::deBoor(int i, int p, double xi, vector<double> KV)

Function description:
This function evaluates a NURBS curve at the given knot locations and returns the
cooresponding coodinates in 2D physical space

Precondition:
A pointer to the structure Geom_data had to be created and passed in


Postcondition:
This function updates the node field within the Geom_data with the physical
coordinate of each xi evaluation point

**********************************************************************************/


double nurb::deBoor(int i, int p, double xi, vector<double> KV)
{


	return 0.0;
}


/**********************************************************************************
Function prototype:
void nurb::make_family(int p)

Function description:
This function creates the deBoor tree so that it only need be traversed later.
This is all in the aim to reduce computation time since basis function evaluation
is the most called function within isogeometric analysis

Precondition:
int i needs to be passed in.  This is largest basis function index which could be asked for. 
int p is the largest polynomial degree of the tree
vector<double> Xi is the Knot vector for the NURBS curve
vector<double> xi is the vector containing the basis function evaluation points


Postcondition:
This function creates the deBoor_fam data structure and returns is so that it can be
used in later traversal.

**********************************************************************************/

void nurb::make_family(int p, int i, vector<double> Xi, vector<double> xi)
{
	int xi_size = xi.size();
	vector<deBoor_fam *> headList;
	for (int xi_cur = 0; xi_cur < xi_size; xi_cur++) {    // this loop goes through all of the evaluation points.   Each of these are on different trees
		deBoor_fam *head_last = new deBoor_fam;
		for (int p_cur = 0; p_cur <= p; p_cur++) {              // this loop starts at p = 0 and slowly goes up each level in the polynomial pyramid

			bool completed = false;
			int i_cur = 0;

			deBoor_fam *head = new deBoor_fam;      // this is just the head of the current level
			deBoor_fam *current = head;

			while (completed == false) {      // this loop works in the same level of the pyramid, figuring out all of the values there. 
										     // can't be for loop because the width is non contant, ie cause its a pyramid
				if (p_cur != 0) {
					head->daughter = head_last;
					head_last->dad = head;
					head->son = head->daughter->brother;
				}

				// special case if p_cur == 0
				if (p_cur == 0) {
					if (Xi[i_cur] <= xi[xi_cur] && xi[xi_cur] < Xi[i_cur + 1]) {
						current->N = 1.0;
					}
					else {
						current->N = 0.0;
					}
					// set all of the values within the struct that aren't suppose to be NULL
					current->p = 0;
					current->i = i_cur;

				}

				// general case,    this is just the Cox deBoor relation
				else {
					double N_left = (xi[xi_cur] - Xi[i_cur]) / (Xi[i_cur + p_cur] - Xi[i_cur])*current->daughter->N;
					if (isnan(N_left)) {
						N_left = 0.0;
					}
					double N_right = ((Xi[i_cur + p_cur + 1] - xi[xi_cur]) / (Xi[i_cur + p_cur + 1] - Xi[i_cur + 1]))*current->son->N;
					if (isnan(N_right)) {
						N_right = 0.0;
					}
					current->N = N_left + N_right;
					current->daughter = current->son->sister;
					current->daughter->dad = current;
					current->i = i_cur;
					current->p = p_cur;

				}
				// determine if the layer of the pyramid is completed and set the value of that basis function to 1 if xi is 1 as well
				if (i_cur == (p - p_cur) + i) {
					completed = true;
					head_last = head;
					if (xi[xi_cur] == 1.0) {
						current->N = 1.0;
					}

				}
				else {   // need to set up the next sibling
					deBoor_fam *next_one = new deBoor_fam;
					next_one->sister = current;
					next_one->daughter = current->son;
					current->brother = next_one;
					next_one->mom = current->dad;
					if (p_cur != 0) {
						next_one->son = next_one->daughter->brother;
					}
					current = next_one;
					i_cur = i_cur++;
				}
			}
		}
		// now perform a push back on the vector containing the heads
		headList.push_back(head_last);
	}
	// push this onto deBoor_tree
	deBoor_tree.push_back(headList);
}


/**********************************************************************************
Function prototype:
void nurb::curve_len(Geom_data *var)

Function description:
This function computes the curve length between each adjacent xi value as well as the
linear edge length


Precondition:
Geom_data *var needs to be passed in.  This is the pointer to the main data struct


Postcondition:
This function updates the curve length field within the struct

**********************************************************************************/

void nurb::curve_len(Geom_data * var)
{
	vector<vector<int>> index = determine_xi_index(var); // now I know the indexes of var->xi which define the bounds of curves
	for (unsigned int i = 0; i < var->node.size()-1; i++) {    // for the plate and hole, there are 9 xi values so there are 8 curves inbetween them.  that is why is go to size-1
	}

}


Bezier_elem * nurb::extraction1D(Geom_data *var)
{
	// First initialize the bezier element struct
	Bezier_elem *head = new Bezier_elem;
	MatrixXd Operator;
	Operator.resize(var->p_deg + 1, var->p_deg + 1);
	Operator.setIdentity();


	int a = var->p_deg + 1;
	int b = a;
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

	return head;
}


/**********************************************************************************
Function prototype:
vector<vector<int>> nurb::determine_xi_index(Geom_data *var)

Function description:
This function determines the indexes of the xi vector to which each arc is described by


Precondition:
Geom_data *var needs to be passed in.  This is the pointer to the main data struct


Postcondition:
This function returns a 2D array containing the xi indexing information
The first column of this array denotes the starting xi index and the second column
denotes the ending index

Each row is for each of the "smooth" bernstein polynomials
what is meant by smooth is that there are no cusps present within it.
If there is a cusp present within it, it is divided into pieces so there is no cusp

**********************************************************************************/
std::vector<std::vector<int>> nurb::determine_xi_index(Geom_data * var)
{
	vector< vector<int> > index;
	int i = 0;
	bool done = false;
	while (done == false){         // since the knot vector is normalized about 1, it will hit this eventually
		if (var->xi[i] == 1)
			done = true;

		int i_end = i + var->p_deg + 1;
		// now determine if there is a p times repetition on the inside of this span'
		vector<vector<double>> saved;
		for (int j = i + 1; j < i_end; j++) {
			if (var->KV[j] != var->KV[i] && var->KV[j] != var->KV[i + var->p_deg + 1]) {
				if (saved.size() == 0) {     // there is nothing in the array
					vector<double> temp(2, 1);
					temp[0] = var->KV[j];
					saved.push_back(temp);
				}
				else if (var->KV[j] == saved.back()[0]) {    // I am at another repeat so increment the count
					saved.back()[1]++;
				}
				else {       // I am at a new number so add a new row to the array
					vector<double> temp(2, 1);
					temp[0] = var->KV[j];
					saved.push_back(temp);
				}

			}
		}   // end of the saved array creation loop

		vector<int> index_temp(3, 0);
		// find where var->KV[i] exists within the var->xi vector  it is guarenteed to be in there by the way we made the xi vector
		int j = 0;
		while (true) {
			if (var->KV[i] == var->xi[j])
				break;
			else
				j++;
		}
		index_temp[0] = j;
		index_temp[2] = i;     // the third column denotes the global basis function index

		// now identify if there are any cases with cusps I need to deal with
		if (saved.size() != 0) {
			for (unsigned int k = 0; k < saved.size(); k++) {
				if (saved[k][1] >= var->p_deg) {        // there is a cusp I need to deal with
					int j = 0;
					while (true) {
						if (saved[k][0] == var->xi[j])
							break;
						else
							j++;
					}
					index_temp[1] = j;      // this is the ending index within the xi vector
					index.push_back(index_temp);
					// now do the next one
					// now deal with getting the starting index beacuse it is a cusp, it is the same as the ending index of the last one
					index_temp[0] = j;
					index_temp[2] = i;
					//now find the ending index of this next segment   because of the fact that there are p+1 repetative knots... if j = 2 and the corresponding location in KV is 0.25 and the next
					// non repeated number is 0.5, index_temp[1] should be what ever the index in xi is 0.5

					int m = 0;
					while (true) {
						if (saved[k][0] == var->xi[m])
							break;
						else
							m++;
					}
					// m is now the index of the first occurance of the repeated knots
					int temp = static_cast<int>(saved[k][1]);    // typecast to get rid of the warnings
					m = m + temp - 1;
					m += var->p_deg;      // now I have the index of the "0.5" as mentioned above

					// find this number in the xi vector
					int n = 0;
					while (true) {
						if (var->KV[m] == var->xi[n])
							break;
						else
							n++;
					}

					index_temp[1] = n;
					index.push_back(index_temp);

					// I am "hard coding these two interations because there is the before cusp and post cusp segments.
				}
			}
		}
		else {
			// find what it would have been otherwise
			int j = 0;
			while (true) {
				if (var->KV[i_end] == var->xi[j])
					break;
				else
					j++;
			}
			index_temp[1] = j;
			index.push_back(index_temp);
			// initialize the next one
			vector<int> index_temp(2, 0);
		}


		i++;

	}
	return index;
}





