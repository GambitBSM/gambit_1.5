/*
 * Nested virtual class example (for SpecBit)
 * 
 * \author Pat Scott
 * \date 2014-06-26
 * 
 */

#include <iostream>
#include "classes.hpp"

int main()
{

  std::cout << std::endl;
  std::cout << "This is the nested virtual class example." << std::endl; 

  specialised_outer<double> dbl_example;
  specialised_outer<int> int_example;

  //Set all the internal values in both objects to the same values.
  dbl_example.Insides1.set_mA(0.5);
  dbl_example.Insides2.set_mB(0.5);
  int_example.Insides1.set_mA(0.5);
  int_example.Insides2.set_mB(0.5);

  //Print the functions of mA and mB, in each object.
  std::cout << dbl_example.Insides1.get_f_of_mA(1.0) << std::endl;
  std::cout << dbl_example.Insides2.get_f_of_mB(1.0, 1.0) << std::endl;
  std::cout << int_example.Insides1.get_f_of_mA(1.1) << std::endl;
  std::cout << int_example.Insides2.get_f_of_mB(1.1, 1.1) << std::endl;

}
