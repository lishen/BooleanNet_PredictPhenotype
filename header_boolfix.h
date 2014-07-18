/*  boolfix: Boolean network fixed-point analyzer (Version 1) 
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


#include <iostream>
#include <iomanip>
#include <bitset>
#include <vector>
#include <map>
#include <fstream>
#include <stdexcept>
#include <string>

using namespace std;

/*=====================================================================================*/
void error(string s)
{
  throw runtime_error(s);
}

/*===============================================================================*/
void Tokenize(const string& str,vector<string>& tokens,const string& delimiters = "\t")
{
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

/*===============================================================================*/
// Read node interaction and node names from file.
void read_interactions(string fname,  vector<vector<short int> > &m, vector<string> &v_nodes)
{
  m.clear();
  v_nodes.clear();

  string line;
  vector<string> tokens;
  ifstream fin;
  map<string,unsigned int> m_node, node_id;
  
  // Read node interaction file and register node names in hash table.
  fin.open(fname.c_str());
  if (! fin) throw runtime_error("Cannot open file");
  while (getline(fin,line)){
     tokens.clear();
     Tokenize(line,tokens,"\t");
     m_node[tokens[0]] = 1;	// node name hash table.
     m_node[tokens[1]] = 1;
  }
  fin.close();
  fin.clear();

  if (m_node.size() != NODES) 
  {
    cerr << m_node.size() << endl;
    map<string, unsigned int>::iterator i;
    for(i = m_node.begin(); i != m_node.end(); i++)
        cerr << i->first << endl;
    error("wrong number of nodes in network interactions file");
  }
 
  // Assign an integer for each node.
  map<string,unsigned int>::iterator it = m_node.begin();
  unsigned int id = 0;

  unsigned int num_nodes = m_node.size();
  v_nodes = vector<string> (num_nodes);	// node name list.

  while(it != m_node.end()){
    node_id[it->first] = id;	// node name->id hash table.
    v_nodes.at(id) = it->first;
    id++;
    ++it;
  }

  // Create node interaction map and initialize it to zeros.
  m = vector<vector<short int> > (num_nodes);	// node interaction map.
  for (int i = 0; i < num_nodes; i++){
    m[i] = vector<short int> (num_nodes);
    for (int j = 0; j < num_nodes; j++){
      m[i][j] = 0;
    }
  } 

  // Open the same file again and read in the interaction types.
  int e_sign;
  fin.open(fname.c_str());
  if (! fin) throw runtime_error("Cannot open file");
  while (getline(fin,line)){
    tokens.clear();
    Tokenize(line,tokens,"\t");
    if (node_id.count(tokens[0])==0 || node_id.count(tokens[1])==0){
      error("node name does not exist");
    }
    if (tokens[2] == "+") {e_sign = 1;}	// activation.
    else if (tokens[2] == "-") {e_sign = -1;}	// repression.
    m.at(node_id[tokens[1]]).at(node_id[tokens[0]]) = e_sign; // row->child, column->parent.
  }
  fin.close();
  fin.clear();
 
}

/*================================================================================*/

void write_interactions(string fname, const vector<vector<short int> > &m, const vector<string> &v_nodes)
{
  ofstream fout;

  fout.open(fname.c_str()); 
  if (! fout) throw runtime_error("Cannot open file");

  for (int i = 0; i < v_nodes.size(); i++){
    fout << "\t" << v_nodes[i];
  }
  fout << endl;
  for (int i = 0; i < v_nodes.size(); i++){
    fout << v_nodes[i];
    for (int j = 0; j < v_nodes.size(); j++){
	fout << "\t" << m.at(i).at(j);
    }
    fout << endl;
  }
  fout.close();
}

/*===============================================================*/

void read_node_property(string fname, const vector<string> &v_nodes, bitset<NODES> &self_degrad, 
			bitset<NODES> &is_ANDnode, vector<int> &clamp)
{
  // file format: node_name \tab is_self_degrad (y/n) \tab is_ANDnode (y/n) \tab clamp/marker (-1/0/1/2/3)
  // -1 no clamp; 0,1 clamp value; 2 seeded zero; 3 seeded one.
  
  self_degrad.reset();
  is_ANDnode.reset();

  map<string,int> m_nodes;	// map a node name to an integer.
  
  for (int i = 0; i < v_nodes.size(); i++){
    m_nodes[v_nodes[i]] = i;	// this should map each node name to the same integer as previous.
  }

  string line;
  vector<string> tokens;
  ifstream fin;
  
  map<string,int> m_test;

  fin.open(fname.c_str());
  if (! fin) throw runtime_error("Cannot open file");
  while (getline(fin,line)){
    tokens.clear();
    Tokenize(line,tokens,"\t");	// split a string into tokens.
    if (! m_nodes.count(tokens[0])) {error("node does not exist");} // match the read-in node with node map already created.
    m_test[tokens[0]] = 1;	// set a tag for each node read in.
	if(m_nodes.find(tokens[0]) != m_nodes.end())
	{
		if (tokens[1] == "y") {self_degrad.set(m_nodes[tokens[0]]);}
		if (tokens[2] == "y") {is_ANDnode.set(m_nodes[tokens[0]]);}
		clamp.at(m_nodes[tokens[0]]) = atoi(tokens[3].data());	// Read clamp/mask values.
	}
  }
  fin.close();
  fin.clear();

  if (m_test.size() != NODES) {error("wrong number of nodes in node property file");}
}

/*==============================================================*/
// Evolve the Boolean network one step forward by sweeping through all nodes.
typedef bitset<NODES> cellstate;

void evolve(const vector<vector<short int> > &m, const bitset<NODES> &self_degrad, 
	const bitset<NODES> &is_ANDnode, const vector<int> &clamp, 
	const cellstate &S_in, cellstate &S_out)
{
	S_out.reset();	// output states are initialized to zeros.

	int mysum;
	int myprod;

	for (int i = 0; i < NODES; i++)
	{
		if(clamp[i] == 0 || clamp[i] == 1)	// if the node is clamped, set it directly to clamped value then skip.
		{
			if(clamp[i] == 1)	// Default = 0.
				S_out.set(i);
			continue;
		}
		if (is_ANDnode.test(i))	// AND node.
		{	
			myprod = 1;
			for (int j = 0; j < NODES; j++)
			{
				if (m.at(i).at(j) == 1 && !S_in.test(j) ||
				m.at(i).at(j) == -1 && S_in.test(j)) // one input zero, output becomes zero.
				{
					myprod = 0;
					break;
				}
			}
			if (myprod == 1)
				S_out.set(i);
			else if(myprod == 0 && !self_degrad.test(i)) 
				S_out[i] = S_in[i];

		}
		else // Normal node.
		{	
			mysum = 0;
			for (int j = 0; j < NODES; j++)
			{
				if(S_in.test(j))
					mysum += m[i][j];
			}
			if (mysum > 0)
				S_out.set(i);       
			else if(mysum == 0 && !self_degrad.test(i)) 
				S_out[i] = S_in[i];
		}
	}
}

/* ==================================================== */
// Check whether state consistent with the clamped values and marker.
bool chkclamp(unsigned long s, const vector<int> &clamp)
{
	cellstate S = cellstate(s);
	for(int i = 0; i < NODES; i++)
	{
		if(clamp[i] == 1 && !S.test(i))	// Clamped to 1 but seed is 0.
			return false;
		else if(clamp[i] == 0 && S.test(i))	// Clamped to 0 but seed is 1.
			return false;
		else if(clamp[i] == 2 && S.test(i))	// seeded zero but seed is 1.
			return false;
		else if(clamp[i] == 3 && !S.test(i))	// seeded one but seed is 0.
			return false;
	}

	return true;
}


/*==================================*/
// Calculate basin counts by evolving Boolean network.
void basin_count(const vector<vector<short int> > &m, const bitset<NODES> &self_degrad, 
		const bitset<NODES> &is_ANDnode, const vector<int> &clamp, 
		map<unsigned long, unsigned long> &basin_coun)
{

  basin_coun.clear();
  cellstate S1, S2;
  S1.set();	// set all bits to one.
  unsigned long max_state = S1.to_ulong();	// number of initializing seeds.
  
  for (unsigned long i = 0; i <= max_state; ++i){
	  if(!chkclamp(i, clamp))	// skip the seed if NOT consistent with clamped values. Therefore, all nodes will start with clamped values properly set.
		  continue;
	  S1 = cellstate(i);	// set initial state.
	  evolve(m,self_degrad,is_ANDnode,clamp,S1,S2);	// S1 evolves to S2.

	  while( S2 != S1 ){	// check convergence.
		  S1 = S2;
		  evolve(m,self_degrad,is_ANDnode,clamp,S1,S2);
	  }
	  basin_coun[S2.to_ulong()]++;
  }  

}

/*==================================*/
// Judge whether the cell state is in sporulation.
bool ispor(cellstate S)
{
	if(S.test(8) && S.test(10))
	//   MMG           Ndt80
		return true;
	return false;
}

/*==================================*/
// Record seeds that lead to sporulation.
void seed_reco(const vector<vector<short int> > &m, const bitset<NODES> &self_degrad, 
	 const bitset<NODES> &is_ANDnode, const vector<int> &clamp, vector<unsigned long> &seed_spor)
{
  cellstate S1, S2;
  S1.set();	// set all bits to one.
  unsigned long max_state = S1.to_ulong();	// number of initializing seeds.
  
  for (unsigned long i = 0; i <= max_state; ++i){
	  if(!chkclamp(i, clamp))	// skip the seed if NOT consistent with clamped values. Therefore, all nodes will start with clamped values properly set.
		  continue;
	  S1 = cellstate(i);	// set initial state.
	  evolve(m,self_degrad,is_ANDnode,clamp,S1,S2);	// S1 evolves to S2.

	  while( S2 != S1 ){	// check convergence.
		  S1 = S2;
		  evolve(m,self_degrad,is_ANDnode,clamp,S1,S2);
	  }
	  if(ispor(S2))
		  seed_spor.push_back(i);
  }  

}

/*=========================================*/

struct gt
{
  bool operator()(unsigned long s1, unsigned long s2) const
  {
    return s1 > s2;
  }
};

/* ======================================== */
// Write counts/percentage of end status to file.
void write_basin_bitset_coun(string fname, map<unsigned long, unsigned long> &basin_coun,
			const vector<string> &v_nodes, const int nclamp = 0)
{
  ofstream fout;
  map<unsigned long, unsigned long>::iterator it = basin_coun.begin();
  multimap<unsigned long, unsigned long,gt> basin_coun_mmap;
  
  while(it != basin_coun.end()){
    basin_coun_mmap.insert(pair<unsigned long, unsigned long>(it->second, it->first));
    it++;
  }

  cellstate S1;
  S1.set();
  unsigned long num_state = S1.to_ulong() + 1;	// total number of initial states.
  num_state >>= nclamp;	// Adjust for the clamped nodes.

  

  fout.open(fname.c_str());
  if (! fout) throw runtime_error("Cannot open file");
  fout << "FixedPoint_ID\t";  
  for (int i = NODES-1; i >=0 ; --i){
    fout << v_nodes[i] << "\t";
  }  
  fout << "basin_size\tbasin_size_percent" << endl;

    

  cellstate S;
  
  multimap<unsigned long, unsigned long,gt>::iterator it_mmap = basin_coun_mmap.begin();

  int id = 1;
  while(it_mmap != basin_coun_mmap.end()){
    S = cellstate(it_mmap->second);
    fout << id << "\t";
    for (int i = NODES-1; i >=0 ; --i){
      fout << S[i] << "\t";
    }
    fout << it_mmap->first << "\t";
    fout << fixed;
    fout << setprecision(4) << 100.0*double(it_mmap->first)/double(num_state) << endl;
    it_mmap++;
    id++;
  }
  fout.close();
}

/* ================================================== */
// Write the seeds that lead to sporulation to file.
void write_seed_spor(string fname, const vector<unsigned long> &seed_spor, 
		const vector<string> &v_nodes, const int nclamp = 0)
{

  ofstream fout;
  
  cellstate S1;
  S1.set();
  unsigned long num_state = S1.to_ulong() + 1;	// total number of initial states.
  num_state >>= nclamp;	// Adjust for the clamped nodes.

  fout.open(fname.c_str());
  if (! fout) throw runtime_error("Cannot open file");
  fout << "Sporulation = " << (double)seed_spor.size()/num_state << endl;
  for (int i = NODES-1; i >0 ; --i){
    fout << v_nodes[i] << "\t";	// nodes are output at reverse order.
  }  
  fout << v_nodes[0] << endl;

  vector<unsigned long> numstat(NODES);	// number statistics for each node.
  for(size_t i = 0; i < seed_spor.size(); i++)	// output all seeds.
  {
	  cellstate S(seed_spor[i]);
	  for (int j = NODES-1; j >0 ; --j)
	  {
		  fout << S[j] << "\t";
		  if(S.test(j))
			  numstat[j]++;
	  }
	  fout << S[0] << endl;
	  if(S.test(0))
		  numstat[0]++;
  }

  // Output number statistics for seeds.
  fout << endl << "<--------- Number statistics for all nodes --------->" << endl;
  for(int i = 0; i < NODES; i++)
	  fout << v_nodes[i] << "\t" << numstat[i] << "\t" << (double)numstat[i]/seed_spor.size() << endl;

}

/*====================================*/
// Number of clamped nodes and markers.
int num_clamped(vector<int> &clamp)
{
	int c = 0;
	for(int i = 0; i < clamp.size(); i++)
	{
		if(clamp[i] != -1)
			c++;
	}
	
	return c;
}

