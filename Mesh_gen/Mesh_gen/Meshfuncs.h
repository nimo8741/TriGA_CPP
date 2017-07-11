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
	std::vector<int> controlP;    // this is a vector which contains the indexes of the nodes that triangular element contains.  It is essentially a row from the IEN array
	std::vector<int> global_side;    // this contains the global index for each of the three sides of the triangle.  They are not unique to each element as almost all sides are shared
	                                 // this index corresponds to the index in the global_edges variable
};

struct quadInfo {
	double quadP[16][2] =      // this contains the x and y coordinates for the quadrature points   this only supports 16 point quadrature rule right now
	// define the values within quadP
	{
		{0.0571, 0.0655},
		{0.2768, 0.0502},
		{0.5836, 0.0289},
		{0.8602, 0.0097},
		{0.0571, 0.3112},
		{0.2768, 0.2386},
		{0.5836, 0.1374},
		{0.8602, 0.0461},
		{0.0571, 0.6317},
		{0.2768, 0.4845},
		{0.5836, 0.2790},
		{0.8602, 0.0936},
		{0.0571, 0.8774},
		{0.2768, 0.6729},
		{0.5836, 0.3875},
		{0.8602, 0.1301}
	};

	double weights[16] = { 0.0236, 0.0354, 0.0226, 0.0054, 0.0442, 0.0663, 0.0423, 0.0102, 0.0442, 0.0663, 0.0423, 0.0102, 0.0236, 0.0354, 0.0226, 0.0054 };  
	                                      // this contains the weighting information
};

struct tri_10_output {
	std::vector<double> R;                    // this is the evaluation of the basis function evaluated at the desired locations
	Eigen::MatrixXd dR_dx;   // this is the array containing the [dR_dx, dR_dy] derivatives for each desired location
	std::vector<double> x;                    // this is the physical location corresponding to the desired parametric location
	double J_det;                             // this is the determinant of the Jacobian which maps from the parametric space to the physical space
	Eigen::Matrix2d Jacob;   // this is the array containing the shape function derivatives which make up the Jacobian
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
		void readMsh(std::string filename, int degree);

		// Part 3: Smoothing and degree elevating the linear mesh
		void smoothMesh(int mesh_degree);
		//void smooth_weights(int degree);
		void create_side_nodes();
		void organize_boundary();
		bool operator () (int i, int j);
		void split_and_extract();
		std::vector<double> determine_elem_split(int cur_nurb, int cur_elem);
		void curve_refine(Bezier_handle * Bez, int cur_elem, std::vector<double> xi_to_add, bool part1);
		void elevate_degree(Bezier_handle * Bez, int element);
		void assign_boundary_points();
		//void boundary_weights(int degree);
		//void refine_and_elevate(int degree, std::vector<double> xi_list, );

		//void adjust_boundary_deg3();
		//std::vector<double> eval_Bez_elem(double xi_val, unsigned int element, unsigned int cur_nurb);
		//void LE2D(int degree);
		//void evaluate_tri_basis(int degree);
		//tri_10_output tri_10_fast(unsigned int elem, int q);

		//Eigen::MatrixXd create_matrix(std::vector<std::vector<double>> input);


		// Part 4: Write the mesh to a .xmsh file to be read in by matlab and displayed
		void create_xmsh(std::string filename, int degree);
		void display_mesh(std::string filename);

	protected:
		std::vector< std::vector<int> > BC;       // this is the boundary condition matrix.  There is only one for the entire problem

		// the length of the vector is the same as the number of NURBS curves.  The size of each cell is the number of xi evaluation points
		std::vector< Bezier_handle *> Elem_list;    // this is a list of he bezier elements.  Each entry contains all of the elements which make up the NURBS curve
		int num_curves;
		std::vector < double> fast_fact;
		std::vector <std::vector< std::vector<int>>> Master_IEN;   // this is the IEN Array for the NURBS curves, not for the triangular elements


		std::vector <Tri_elem *> triangles;      // this is a list of all of the triangular elements in the mesh
		std::vector<std::vector <int>> global_edges; 
		// the first column of the above field denotes the index of the first node in the side, 
		// the second is the second node for the side, and the third is the physical group
	    // to which that edge belongs to.  By physical side, it is meant which bezier "half element" it belongs to
		std::vector<std::vector <double>> node_list;
		std::vector<int> bNodes;    // this is a list of all of the boundary nodes
		int phy_groups;     // This is a value which keeps track of how many boundary line segments there are.  It is 1 larger than that since the last physical group is the internal domain of interest
		std::vector<std::vector<double>> tri_N;                 // this is the array containing the parametric evaluations.  The first dimension is each of the 16 evaluation points and the second is each of the 10 basis functions
	    std::vector<std::vector<std::vector<double>>> tri_dN_du;
		// this is the 3D array containing the parametric derivative information.  The first dimension is the 16 evaluation points, the second are the 10 basis functions, and the third is each of the 3 directions of derivatives.
		
		int degree;
		int nodes_in_triangle;
		std::vector<std::vector<unsigned int>> tri_NURB_elem_section_side;     // this is a 2D array while holds the information for the triangles along the boundary, the NURBS curve of the boundary, the element of this curve, and the section of the element
		std::vector<std::vector<int>> node_side_index;




	private:
		void refine_Xi(Geom_data *var);
		const int num_quad = 16;



};




#endif // MESHFUNCS_H