/*
    boolfix: Boolean network fixed-point analyzer (version 1)
    Copyright (C) 2008  Iouri Chepelev

	Modified by: Li Shen
	Date: August, 2008

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <time.h>
#include "header_boolfix.h"

int main(int argc, char* argv[]){

  if (argc < 4) error("Input format must be: interaction_filename  node_property_filename output_name");
  
  string in_fname = string(argv[1]); 
  string node_prop_fname = string(argv[2]); 
  string matrix_fname = string(argv[3]) + "_matrix.txt";
  string fixedpoint_fname = string(argv[3]) + "_fixedpoints.txt";

  vector<vector<short int> > m;	// node interaction map.
  vector<string> v_nodes;	// node names.
  read_interactions(in_fname,m,v_nodes);
  write_interactions(matrix_fname,m,v_nodes);

  bitset<NODES> self_degrad;
  bitset<NODES> is_ANDnode;
  vector<int> clamp(NODES, -1);	// set clamp values for nodes. default= -1 (NO clamp).
  read_node_property(node_prop_fname,v_nodes,self_degrad,is_ANDnode,clamp);
 
  // Boolean network kicks off.
  time_t start, end;
  time(&start);
  map<unsigned long, unsigned long> basin_coun;  
  basin_count(m,self_degrad,is_ANDnode,clamp,basin_coun);
  write_basin_bitset_coun(fixedpoint_fname,basin_coun,v_nodes,num_clamped(clamp));
  time(&end);
  double secs = difftime(end, start);
  cerr << "Boolean network cost: " << secs << " seconds." << endl; 
 
  return 0;
}
