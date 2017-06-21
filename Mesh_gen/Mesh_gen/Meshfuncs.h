#ifndef MESHFUNCS_H
#define MESHFUNCS_H
#include <vector>
#include <Eigen/Dense>

struct Geom_data {
	std::vector< std::vector<double> > controlP;   // The reason there are triple points are so that the first two dimensions are X and Y points  
	std::vector< int > B_flag;     // Similar to controlP, the first is the side number and the second is the boundary flag ID   
	std::vector<double> KV;          // This is a 1D array containing the data for the knot vector  
	std::vector<double> weight;     // This is a 1D array containing the data for the weights of each control point
	std::vector< std::vector<double> > node;      // This is the 2D array which contains the physical xy location of each of the xi evaluation points
	Geom_data *next;   // I'm trying to add some sort of linked list to be able to traverse between all of the distinct NURBS curve elements.  Might want to go with something else entirely
	Geom_data *previous;
	int p_deg;     // polynomial degree of the spline
	int NCP;     // number of control points
	int LKV;     // length of the spline's knot vector
	
};



struct Bezier_elem {    // This is the list structure that will eliminate the need to create the tree for every evaluation.
	std::vector <std::vector<double>> controlP; // this array holds the data for the control points and weights.  the x location is in the first column, y in the second, and weight in the third (for when projected)
	std::vector<double> weight;    // this holds the weight of the post-extraction local Bezier control points
	std::vector<std::vector<double>> xi_real_loc;        // this is the physical location of where xi is evaluated at
	double edge_len = 0.0;                  // this is the linear edge length along the entire NURBS curve
	double curve_len = 0.0;                 // this is the approximate length of the actual curve length
	int side;                           // this is the side of the global NURBS curve in which this element lives

};

struct Bezier_handle {
	std::vector <Eigen::MatrixXd> Operator;
	int p;
	int n_el;
	std::vector <std::vector <double>> N;                 // this is the value of the evaluate basis function each row is a different evaluation point and each col is a new basis function index
	std::vector <double> xi_evals;           // this is the vector of xi evaluation points
	std::vector<Bezier_elem> elem_Geom;     // this is a vector which contains all of the data for the control points and weights
};

struct Tri_elem {
	std::vector<std::vector<double>> controlP;    // this is a 2D array containing the x and y location of each of the 10 control points.
	std::vector<int> global_side;    // this contains the global index for each of the three sides of the triangle.  They are not unique to each element as almost all sides are shared
	                                 // this index corresponds to the index in the global_edges variable
	
};



class nurb{
	public:
		nurb();
		//~nurb();   // this is the deconstructer I think but I'm not sure if I'll need it

		// Part 1: reading in, extracting, and refining the NURBS curves
		Geom_data* readSpline(std::string filename);
		void rearrange(Geom_data *var, int array, int B_rows);
		bool hasnumber(std::string line);
		void NURBS2poly(Geom_data *var);
		void evalBez(Geom_data *var, Bezier_handle * Bez, int elem, std::vector<double> KV_old);
		Geom_data *head_nurb;                               // this is the head of the linked list of the structs that contain the nurb data
		void curve_len(Bezier_handle *Bez, int i);
		Bezier_handle * extraction1D(Geom_data *var);
		void fact_table();
		void eval_bern(int p, int n, Bezier_handle *Bez);
		void get_P_and_W(Bezier_handle *Bez, Geom_data *var, std::vector<double> KV_old);
		std::vector<std::vector<int>> IEN(int n, int p, std::vector<double> Xi,int n_el);
		void subdivide_element(Bezier_handle *Bez, int e_index);



		// Part 2: Making the linear mesh and reading the result back into this program
		void create_geo_file(std::string filename);
		void call_gmsh(std::string filename);
		void readMsh(std::string filename);

		// Part 3: Smoothing and degree elevating the linear mesh
		void smoothMesh();


	protected:
		std::vector< std::vector<int> > BC;       // this is the boundary condition matrix.  There is only one for the entire problem

		// the length of the vector is the same as the number of NURBS curves.  The size of each cell is the number of xi evaluation points
		std::vector< Bezier_handle *> Elem_list;    // this is a list of he bezier elements.  Each entry contains all of the elements which make up the NURBS curve
		int num_curves;
		std::vector < double> fast_fact;
		std::vector <std::vector< std::vector<int>>> Master_IEN;


		std::vector <Tri_elem *> triangles;      // this is a list of all of the triangular elements in the mesh
		std::vector<std::vector <int>> global_edges; 
		// the first column of the above field denotes the index of the first node in the side, 
		// the second is the second node for the side, and the third is the physical group
	    // to which that edge belongs to.  By physical side, it is meant which boundary condition group it belongs to.
		std::vector<std::vector <double>> node_list;





	private:
		void refine_Xi(Geom_data *var);


};




#endif // MESHFUNCS_H