#ifndef MESHFUNCS_H
#define MESHFUNCS_H
#include <vector>
#include <Eigen/Dense>

struct Geom_data {
	std::vector< std::vector<double> > controlP;   // The reason there are triple points are so that the first two dimensions are X and Y points  
	std::vector< std::vector<int> > B_flag;     // Similar to controlP, the first is the side number and the second is the boundary flag ID   
	std::vector<double> KV;          // This is a 1D array containing the data for the knot vector  
	std::vector<double> weight;     // This is a 1D array containing the data for the weights of each control point
	std::vector<double> xi;         // This is the 1D array of xi values at which to evaluate the NURBS curves
	std::vector< std::vector<double> > node;      // This is the 2D array which contains the physical xy location of each of the xi evaluation points
	std::vector<double> edge_length;      // each entry in this vector is the total edge length by linearally putting a line between each consecutive evaluation point
	std::vector<double> curve_length;      // same as above but of the analytical curve length along the NURBS curve
	Geom_data *next;   // I'm trying to add some sort of linked list to be able to traverse between all of the distinct NURBS curve elements.  Might want to go with something else entirely
	Geom_data *previous;
	int p_deg;     // polynomial degree of the spline
	int NCP;     // number of control points
	int LKV;     // length of the spline's knot vector
	
};



/*
//////////////////////////////////////////////////////////////////////
// Need to add a Bezier Element object and then do 
extraction on the original NURBS (in plate and hole, I will get two of
these objects).  From here I then can do all of the stuff with the 
bezier objects.  The deBoor tree isnt really needed since all of those
quadrature point values can be found using the simple N choose K notation
and saved in a linked list.  These will be the same values for each of the elements
The only thing I will have to be careful with is the data containing the 
control points and the weights
**********************************************************************
*/
struct deBoor_fam {    // This is the tree structure that will eliminate the need to create the tree for every evaluation.  Simply create it once and then just traverse
	int i;                    // this is the index of the basis function
	int p;                    // this is the polynomial degree
	double N;                 // this is the value of the evaluate basis function

	/*     Girls are to the left, guys are to the right             */
	
	deBoor_fam *daughter = NULL;     // this is equivalent to saying the left child
	deBoor_fam *son = NULL;          // this is equivalent to saying the right child

	deBoor_fam *brother = NULL;      // this is equivalent to saying the right sibling
	deBoor_fam *sister = NULL;       // this is equivalent to saying the left sibling

	deBoor_fam *mom = NULL;          // this is equivalent to saying the left parent
	deBoor_fam *dad = NULL;          // this is equivalent to saying the right parent
	
};

struct Bezier_elem {
	std::vector <Eigen::MatrixXd> Operator;

	int n_el;


};

class nurb{
	public:
		nurb();
		//~nurb();   // this is the deconstructer I think but I'm not sure if I'll need it
		Geom_data* readSpline(std::string filename);
		void rearrange(Geom_data *var, int array, int B_rows);
		bool hasnumber(std::string line);
		void NURBS2poly(Geom_data *var);
		void evalNURBS(Geom_data *var);
		double deBoor(int i, int p, double xi, std::vector<double> KV);
		void make_family(int p, int i, std::vector<double> Xi, std::vector<double> xi);
		Geom_data *head_nurb;                               // this is the head of the linked list of the structs that contain the nurb data
		void curve_len(Geom_data *var);
		Bezier_elem * extraction1D(Geom_data *var);



	protected:
		std::vector< std::vector<int> > BC;       // this is the boundary condition matrix.  There is only one for the entire problem
		std::vector< std::vector<deBoor_fam *> > deBoor_tree;           // this is a vector which will contain the heads of each of the deBoor families.
		// the length of the vector is the same as the number of NURBS curves.  The size of each cell is the number of xi evaluation points
		std::vector< Bezier_elem *> Elem_list;    // this is a list of the heads of the bezier elements.  Each head starts the list for its respective NURBS curve.
		int num_curves;


	private:
		std::vector< std::vector<int> > determine_xi_index(Geom_data *var);


};




#endif // MESHFUNCS_H