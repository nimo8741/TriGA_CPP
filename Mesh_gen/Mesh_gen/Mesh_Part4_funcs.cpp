#include "Meshfuncs.h"
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <cfloat>
#include <time.h>
#include <engine.h>

#pragma comment ( lib, "libmat.lib" )
#pragma comment ( lib, "libmx.lib" )
#pragma comment ( lib, "libmex.lib" )
#pragma comment ( lib, "libeng.lib" )

using namespace std;

void nurb::create_xmsh(string filename, int degree) 
{
	ofstream xmsh_file;
	xmsh_file.open(filename + ".xmsh");
	
	// get system time
	time_t cur_time = time(NULL);
	char str[26];
	ctime_s(str, sizeof str, &cur_time);


	xmsh_file << "% MATLAB output file for " << filename << endl;
	xmsh_file << "% Completed at " << str << endl;
	xmsh_file << "% Degree of mesh:" << endl << degree << endl << endl;

	// loop through the node list

	xmsh_file << "% Node List" << endl << "% NUMBER OF NODES" << endl << node_list.size() << endl;
	unsigned int node_size = node_list.size();
	for (unsigned int i = 0; i < node_size; i++) {
		xmsh_file << i+1 << "," << node_list[i][0] << "," << node_list[i][1] << "," << node_list[i][2] << endl;
	}
	xmsh_file << "% END OF NODE LIST " << endl << endl;

	// loop through the connectivity array

	xmsh_file << "% Connectivity array" << endl << "% NUMBER OF TRIANGLES, NODES PER TRIANGLE" << endl << triangles.size() << "," << nodes_in_triangle << endl;
	unsigned int tri_size = triangles.size();
	for (unsigned int i = 0; i < tri_size; i++) {
		xmsh_file << i + 1 << ",";
		for (int j = 0; j < nodes_in_triangle; j++) {
			if (j == nodes_in_triangle - 1)
				xmsh_file << triangles[i]->controlP[j] + 1 << endl;
			else
				xmsh_file << triangles[i]->controlP[j] + 1 << ",";
		}
	}
	xmsh_file << "% END OF CONNECTIVITY ARRAY" << endl;
}

void nurb::display_mesh(string filename)
{
	//Engine *mat_eng;
	//mat_eng = engOpen("null");
}