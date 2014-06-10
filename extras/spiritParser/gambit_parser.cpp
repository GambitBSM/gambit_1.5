//////////////////////////////////////////////////////////
// Simple example for GAMBIT ini-file parser
// 
// Christoph Weniger, 6 May 2013
//////////////////////////////////////////////////////////

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
// #include <boost/config/warning_disable.hpp>
// #include <boost/spirit/include/phoenix_core.hpp>
// #include <boost/spirit/include/phoenix_operator.hpp>
// #include <boost/spirit/include/phoenix_object.hpp>

#include <iostream>
#include <fstream>
#include <string>

// Define structures to store ini-file data
namespace client
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    struct backend
    {
      std::string name;
      std::string backend;
    };

    struct dependency
    {
      std::string name;
      std::string module;
    };

    struct observable
    {
      std::string name;
      std::string module;
      std::vector<backend> backends;
      std::vector<dependency> dependencies;
    };
}

// Tell boost fusion how these stuctures look
BOOST_FUSION_ADAPT_STRUCT(
    client::observable,
    (std::string, name)
    (std::string, module)
    (std::vector<client::backend>, backends)
    (std::vector<client::dependency>, dependencies)
)

BOOST_FUSION_ADAPT_STRUCT(
    client::dependency,
    (std::string, name)
    (std::string, module)
)

BOOST_FUSION_ADAPT_STRUCT(
    client::backend,
    (std::string, name)
    (std::string, backend)
)

namespace client
{
    // Define grammar for backend
    template <typename Iterator>
    struct backend_parser : qi::grammar<Iterator, backend()>
    {
        backend_parser() : backend_parser::base_type(start)
        {
            using qi::lit;
            using ascii::char_;
            using boost::spirit::eol;

            space %= char_(" ");
            valid_string %= char_("a-zA-Z_") >> *(char_("a-zA-Z0-9_."));
            keyword %= lit("-Backend");

            start %=
                *space
                >> keyword 
                >> +space
                >> valid_string >> +space 
                >> valid_string >> *space 
                >> eol
                ;
        }

        qi::rule<Iterator, std::string()> valid_string;
        qi::rule<Iterator> space;
        qi::rule<Iterator> keyword;
        qi::rule<Iterator, backend()> start;
    };

    // Define grammar for dependency
    template <typename Iterator>
    struct dependency_parser : qi::grammar<Iterator, dependency()>
    {
        dependency_parser() : dependency_parser::base_type(start)
        {
            using qi::lit;
            using boost::spirit::eol;
            using ascii::char_;

            space %= char_(" ");
            valid_string %= char_("a-zA-Z_") >> *(char_("a-zA-Z0-9_."));
            keyword %= lit("-Dependency");

            start %=
                *space
                >> keyword 
                >> +space
                >> valid_string >> +space 
                >> valid_string >> *space 
                >> eol
                ;
        }

        qi::rule<Iterator, std::string()> valid_string;
        qi::rule<Iterator> space;
        qi::rule<Iterator> keyword;
        qi::rule<Iterator, dependency()> start;
    };

    // Define grammar for observable entry (including backend/dependency)
    template <typename Iterator>
    struct observable_parser : qi::grammar<Iterator, observable()>
    {
        observable_parser() : observable_parser::base_type(start)
        {
            using qi::lit;
            using boost::spirit::eol;
            using ascii::char_;


            space %= char_(" ");
            valid_string %= char_("a-zA-Z_") >> *(char_("a-zA-Z0-9_."));
            keyword %= lit("Observable");

            start %=
                *space
                >> keyword >> +space
                >> valid_string >> +space
                >> valid_string >> *space
                >> eol
                >> *my_backend
                >> *my_dependency
                ;
        }

        client::dependency_parser<Iterator> my_dependency;
        client::backend_parser<Iterator> my_backend;
        qi::rule<Iterator, std::string()> valid_string;
        qi::rule<Iterator> space;
        qi::rule<Iterator> keyword;
        qi::rule<Iterator, observable()> start;
    };
}

int main()
{
    typedef std::string::const_iterator iterator_type;
    typedef client::observable_parser<iterator_type> observable_parser;

    observable_parser g; // Our main parser
    std::vector<client::observable> result;

    // Read in ini file
    std::ifstream in("gambit.ini");
    std::string str((std::istreambuf_iterator<char>(in)), 
        std::istreambuf_iterator<char>());

    // Parse it!
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bool r = parse(iter, end, *g, result);

    // And pretty print what is now in the structs
    for(std::vector<client::observable>::iterator it = result.begin(); it !=
        result.end(); ++it)
    {
      std::cout << "Observable (Module): " << (*it).name << " (" << (*it).module << ")" << std::endl;
      for(std::vector<client::backend>::iterator it2 = (*it).backends.begin(); it2 !=
          (*it).backends.end(); ++it2)
      {
        std::cout << "-Backend (Module): " << (*it2).name << " (" << (*it2).backend << ")" << std::endl;
      }
      for(std::vector<client::dependency>::iterator it2 = (*it).dependencies.begin(); it2 !=
          (*it).dependencies.end(); ++it2)
      {
        std::cout << "-Dependency (Module): " << (*it2).name << " (" << (*it2).module << ")" << std::endl;
      }
    }

    // done!
    return 0;
}
