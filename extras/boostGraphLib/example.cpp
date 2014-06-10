//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
//
//  Towards dependency resolution with the
//  boost graph library
//
//  *********************************************
//
//  Authors
//  =======
//
//  (add name and date if you modify)
//
//  Christoph Weniger
//  Jan 29 2012
//
//  *********************************************

#include <boost/assign/std/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <queue>
#include <fstream>
#include <boost/graph/graphviz.hpp>

using namespace boost::assign;
using namespace boost;
using namespace std;

// Information associated with vertices and edges of boost graph:
struct Vertex {
  Vertex()
  {
    function = "nothing";
  }
  vector<string> in; // list of input variables
  vector<string> out; // list of output variables
  string module; // module name
  string function; // function name
  int status; // function disabled: 0 , available: 1, active: 2
  // TODO: 
  // - Add pointer to associated function
  // - Add method to pass correct information to/from masterDict
};
struct Edge {
  string variable; // name of associated variable
};

// Central graph object and its types
typedef adjacency_list<vecS, vecS, directedS, Vertex, Edge> MasterGraphType;
typedef graph_traits<MasterGraphType>::vertex_descriptor VertexID;
typedef graph_traits<MasterGraphType>::edge_descriptor EdgeID;
MasterGraphType masterGraph;

// Convenience macros to construct vertices
#define NEW_NODE(MODULE, FUNCTION) \
  current_vertex = add_vertex(masterGraph); \
  masterGraph[current_vertex].module = MODULE; \
  masterGraph[current_vertex].function = FUNCTION; \
  masterGraph[current_vertex].status = 1;
#define DISABLE_NODE() \
  masterGraph[current_vertex].status = 0;
#define ADD_INPUT(...) \
  masterGraph[current_vertex].in += __VA_ARGS__;
#define ADD_OUTPUT(...) \
  masterGraph[current_vertex].out += __VA_ARGS__;

// Exemplary vertex initialization (a poor man's meal)
VertexID initialize_vertices() {
  cout << "Initializing vertices: " << endl;
  VertexID current_vertex;
  // TODO: This should happen during rollcall or whatever
  NEW_NODE("Sauce", "roast")
    ADD_INPUT("onions", "heat")
    ADD_OUTPUT("roasted_onions")
  NEW_NODE("Sauce", "lazy_sauce")
    ADD_INPUT("readymeal")
    ADD_OUTPUT("sauce")
    DISABLE_NODE();
  NEW_NODE("Sauce", "mix");
    ADD_INPUT("roasted_onions", "salt", "pepper", "tomatoes");
    ADD_OUTPUT("proto_sauce");
  NEW_NODE("Sauce", "cook");
    ADD_INPUT("proto_sauce", "heat");
    ADD_OUTPUT("sauce");
  NEW_NODE("Pasta", "cook");
    ADD_INPUT("pasta", "water", "heat", "salt");
    ADD_OUTPUT("proto_pasta_aldente", "pasta_water");
  NEW_NODE("Pasta", "remove_water");
    ADD_INPUT("proto_pasta_aldente", "pasta_water");
    ADD_OUTPUT("pasta_aldente");
  NEW_NODE("Dessert", "tiramisu");
    ADD_INPUT("biscuit", "tomatoes", "pepper");
    ADD_OUTPUT("tiramisu");

  // And we need one input node...
  NEW_NODE("CoreBit", "Alpha");
    ADD_OUTPUT("heat", "pasta", "water", "salt", "pepper", "tomatoes", "onions");
  // ...and one output node...
  NEW_NODE("CoreBit", "Omega");
    ADD_INPUT("sauce", "pasta_aldente");
  // ...which specify initial and final variables
  // Omega is required as seed for edge initialization
  return current_vertex; 
}

// Fill wishlist with in-variables from vertex
void fill_wishlist(queue<pair<string, VertexID> > *wishlist, VertexID vertex) {
  for (vector<string>::iterator it = masterGraph[vertex].in.begin(); it !=
      masterGraph[vertex].in.end(); ++it) {
    (*wishlist).push(*(new pair<string, VertexID> (*it, vertex)));
  }
}

// Create map from variable names to the vertices that could actually provide
// them
multimap<string, VertexID> initialize_variableMap() {
  multimap<string, VertexID> variableMap;
  graph_traits<MasterGraphType>::vertex_iterator vi, vi_end;
  for (tie(vi, vi_end) = vertices(masterGraph); vi != vi_end; ++vi) {
    for (vector<string>::iterator si = masterGraph[*vi].out.begin(); si != masterGraph[*vi].out.end(); ++si) {
      if (masterGraph[*vi].status > 0) 
        variableMap.insert(*(new pair<string, VertexID> (*si, *vi)));
    }
  }
  return variableMap;
}

// Print list of vertices
void list_vertices() {
  graph_traits<MasterGraphType>::vertex_iterator vi, vi_end;
  for (tie(vi, vi_end) = vertices(masterGraph); vi != vi_end; ++vi) {
    cout << "  " << masterGraph[*vi].module << "." << masterGraph[*vi].function << endl;
  }
  cout << endl;
}

// Print content of variableMap
void list_variableMap(multimap<string, VertexID> variableMap) {
  cout << "List of available output variables: " << endl;
  for (multimap<string, VertexID>::iterator it = variableMap.begin(); it != variableMap.end(); ++it) {
    cout << "  " << (*it).first << " (" << masterGraph[(*it).second].module << "." ;
    cout << masterGraph[(*it).second].function << ")" << endl;
  }
  cout << endl;
}

// Recursively initialize edges that are required in order to get Omega
void initialize_edges(queue<pair<string, VertexID> > wishlist, multimap<string, VertexID> variableMap) {
  cout << "Initializing edges:" << endl;
  int key_multiplicity;
  bool ok;
  VertexID fromVertex, toVertex;
  string var;
  EdgeID current_edge;
  while (not wishlist.empty()) {
    var = wishlist.front().first;
    toVertex = wishlist.front().second;
    cout << "  " << var << ": ";
    key_multiplicity = variableMap.count(var);
    if ( key_multiplicity == 0 ) {
      cout << "Cannot resolve dependency." << endl;
      exit(0);
    }
    if ( key_multiplicity > 1 ) {
      cout << "!!provided by multiple functions!!, " << endl;
      // TODO: Add some ini-file based resolution 
      // (for now we just take randomly the first matching vertex)
    }
    fromVertex = (*variableMap.find(var)).second;
    if ( masterGraph[fromVertex].status != 2 ) {
      masterGraph[fromVertex].status = 2;
      fill_wishlist(&wishlist, fromVertex);
    }
    cout << masterGraph[fromVertex].module << "." << masterGraph[fromVertex].function << " --> ";
    cout << masterGraph[toVertex].module << "." << masterGraph[toVertex].function << endl;
    tie(current_edge, ok) = add_edge(fromVertex, toVertex, masterGraph);
    masterGraph[current_edge].variable = var;
    wishlist.pop();
  }
  cout << endl;
}

list<int> run_topological_sort() {
  list<int> topo_order;
  topological_sort(masterGraph, front_inserter(topo_order));
  return topo_order;
}

void execute_functions(list<int> topo_order) {
  cout << "Let the games begin (in that order):" << endl;
  for(list<int>::const_iterator i = topo_order.begin();
      i != topo_order.end();
      ++i)
  {
    cout << "  " << masterGraph[*i].module << "." << masterGraph[*i].function;
    if ( masterGraph[*i].status == 0 ) cout << " (disabled)";
    if ( masterGraph[*i].status == 1 ) cout << " (available)";
    if ( masterGraph[*i].status == 2 ) cout << " (ACTIVE)";
    cout << endl;
  }
}

int main() {
  Vertex my_vertex;
  VertexID omega_vertex;
  queue<pair<string, VertexID> > wishlist;
  multimap<string, VertexID> variableMap;
  list<int> topo_order;
  omega_vertex = initialize_vertices();
  omega_vertex = add_vertex(my_vertex, masterGraph);
  std::cout << masterGraph[omega_vertex].function << std::endl;
  return 0;
  list_vertices();
  fill_wishlist(&wishlist, omega_vertex);
  variableMap = initialize_variableMap();
  list_variableMap(variableMap);
  initialize_edges(wishlist, variableMap);
  topo_order = run_topological_sort();
  execute_functions(topo_order);
  ofstream outf("out.gv");
  write_graphviz(outf, masterGraph,
      make_label_writer(get(&Vertex::function, masterGraph)),
      make_label_writer(get(&Edge::variable, masterGraph)));
}
