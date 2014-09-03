// Example for how to use the YAML-cpp parser for the GAMBIT ini file
//
// Christoph Weniger, 2013-06-01
//

#include <fstream>
#include "yaml-cpp/yaml.h"

// Define storage structures
struct dependency
{
  std::string capability;
  std::string module;
};

struct observable
{
  std::string scannerID;
  std::string capability;
  std::string backend;
  std::vector<dependency> dependencies;
};

struct parameter
{
  std::string name;
  std::pair<double, double> range;
  std::string prior;
};

YAML::Node flagNode;

template<typename TYPE> TYPE getFlag(std::string key)
{
  // std::cout << flagNode.size() << std::endl;
  return flagNode[key].as<TYPE>();
}

namespace YAML {
  template<> struct convert<observable>
  {
    static bool decode(const Node& node, observable& rhs)
    {
      rhs.scannerID = node["scannerID"].as<std::string>();
      rhs.capability = node["capability"].as<std::string>();
      rhs.backend = node["backend"].as<std::string>();
      for(YAML::const_iterator it=node["dependencies"].begin();
          it!=node["dependencies"].end(); ++it)
      {
        rhs.dependencies.push_back((*it).as<dependency>());
      }
      return true;
    };
  };
  template<> struct convert<dependency>
  {
    static bool decode(const Node& node, dependency& rhs)
    {
      rhs.capability = node["capability"].as<std::string>();
      rhs.module = node["module"].as<std::string>();
      return true;
    };
  };
  template<> struct convert<parameter>
  {
    static bool decode(const Node& node, parameter& rhs)
    {
      rhs.name =  node["parameter"].as<std::string>();
      rhs.range = node["range"].as<std::pair<double,double> >();
      rhs.prior = node["prior"].as<std::string>();
      return true;
    };
  };
};

void test1()
{
  // Read ini-file file
  std::vector<YAML::Node> roots = YAML::LoadAllFromFile("gambit.yaml");

  // Set central nodes
  YAML::Node inputNode = roots[0];
  YAML::Node outputNode = roots[1];
  flagNode = roots[2];

  // Read observables
  std::vector<observable> observables;
  for(YAML::const_iterator it=outputNode.begin(); it!=outputNode.end(); ++it)
  {
    observables.push_back((*it).as<observable>());
  }

  // Read scanner parameters
  std::vector<parameter> parameters;
  for(YAML::const_iterator it=inputNode.begin(); it!=inputNode.end(); ++it)
  {
    parameters.push_back((*it).as<parameter>());
  }
}


void printRound(YAML::Node node, int depth = 0)
{
    if (node.Type() == YAML::NodeType::Scalar)
    {
        for (int j = 0; j < depth; j++) std::cout << " ";
        std::cout << node << std::endl;
        return;
    }
    if (node.Type() == YAML::NodeType::Sequence)
    {
        for (unsigned int i = 0; i<node.size(); ++i)
        {
            printRound(node[i], depth + 1);
        }
        return;
    }
    if (node.Type() == YAML::NodeType::Map)
    {
        for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
        {
            printRound(it->second, depth + 1);
        }
        return;
    }
}

int importRound(YAML::Node node)
{
    int counter = 0;
    if (node.Type() == YAML::NodeType::Scalar)
    {
        if ( node.Tag() == "!import" )
        {
            std::string filename = node.as<std::string>();
            YAML::Node imported = YAML::LoadFile(filename);
            node = imported;
            std::cout << "Importing: " << filename << std::endl;
            return 1;
        }
    }
    if (node.Type() == YAML::NodeType::Sequence)
    {
        for (unsigned int i = 0; i<node.size(); ++i)
        {
            counter += importRound(node[i]);
        }
        return counter;
    }
    if (node.Type() == YAML::NodeType::Map)
    {
        for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
        {
            counter += importRound(it->second);  // Only values are processed
        }
        return counter;
    }
}

void test_import()
{
    YAML::Node node = YAML::LoadFile("master.yml");
    printRound(node);
    int import_counter = 0;
    int last_import_counter = 0;
    for ( int i = 0; i < 10; ++i)
    {
        last_import_counter = importRound(node);
        import_counter += last_import_counter;
    }
    if (last_import_counter > 0)
    {
        std::cout << "Recursion depth for YAML imports exceeded." << std::endl;
        exit(1);
    }
    printRound(node);
    std::cout << "Number of YAML imports: " << import_counter << std::endl;
}

int main()
{
    //test1();
    test_import();
    return 0;
}
