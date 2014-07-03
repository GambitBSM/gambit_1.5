#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char **argv)
{
        if (argc == 1)
        {
                std::cout << "fortran72:  No input files." << std::endl;
                return 0;
        }
        else if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "-help" || std::string(argv[1]) == "--h" || std::string(argv[1]) == "--help")
        {
                std::cout << "This program ensures that fortran 77 files stay with the 72 column limit.\n"
                        << "Usage:  fortran72: [files]." << std::endl;
                return 0;
        }
        
        for (int i = 1; i < argc; i++)
        {
                std::ifstream in(argv[i]);
                if (in.is_open())
                {
                        std::cout << "File:  " << argv[i] << std::flush;
                        std::vector<std::string> strs;
                        std::string str;
                        while(getline(in, str))
                        {
                                strs.push_back(str);
                        }
                        
                        std::ofstream out(argv[i]);
                        int line = 0, changed = 0;
                        for (std::vector<std::string>::iterator it = strs.begin(), end = strs.end(); it != end; it++)
                        {
                                line++;
                                size_t pos = (*it + "/0").find_first_of('!');
                                str = it->substr(0, pos);
                                if (str.size() > 72 && str.substr(0, 1) != "*")
                                {
                                        changed++;
                                        std::string comment;
                                        if (pos == std::string::npos)
                                                comment = "";
                                        else
                                                comment = it->substr(pos);
                                        
                                        int lchanged = 0;
                                        while (str.size() > 72)
                                        {
                                                out << str.substr(0, 72) << std::endl;
                                                str = "     &" + str.substr(72);
                                                if (lchanged > 0)
                                                        std::cout << "\n   line " << line << " added." << std::flush;
                                                line++, lchanged++;
                                        }
                                        
                                        if (str.size() > 6)
                                        {
                                                out << str + comment << std::endl;
                                                std::cout << "\n   line " << line << " added." << std::flush;
                                                line++;
                                        }
                                        else
                                        {
                                                out << comment << std::endl;
                                                std::cout << "\n   line " << line << " added." << std::flush;
                                                line++;
                                        }
                                }
                                else
                                {
                                        out << *it << std::endl;
                                }
                        }
                        
                        if (changed == 0)
                                std::cout << " ...OK..." << std::flush;
                        
                        std::cout << std::endl;
                        
                        out.close();
                }
                else
                {
                        std::cout << "File " << argv[i] << " does not exist." << std::endl;
                }
                in.close();
        }
        
        return 0;
}