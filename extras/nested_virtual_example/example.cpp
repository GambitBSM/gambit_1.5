/*
 * Nested virtual class example (for SpecBit)
 * 
 * \author Pat Scott
 * \date 2014-06-26
 * 
 */

#include <iostream>
//#include "classes.hpp"
#include "classes_nonnested.hpp"

int main()
{

  std::cout << std::endl;
  std::cout << "This is the nested virtual class example." << std::endl; 

  specialised_outer<double> dbl_example;
  specialised_outer<int> int_example;

  //Set all the internal values in both objects to the same values.
  dbl_example.set_somevar(2);
  dbl_example.Insides1.set_mA(0.4);
  dbl_example.Insides2.set_mB(0.4);
  int_example.set_somevar(2);
  int_example.Insides1.set_mA(0.4);
  int_example.Insides2.set_mB(0.4);

  //Print the functions of mA and mB, in each object.
  std::cout << dbl_example.Insides1.get_f_of_mA(1.1) << std::endl;
  std::cout << dbl_example.Insides2.get_f_of_mB(1.1, 1.1) << std::endl;
  std::cout << int_example.Insides1.get_f_of_mA(1.1) << std::endl;
  std::cout << int_example.Insides2.get_f_of_mB(1.1, 1.1) << std::endl;

  //Create some references-to-base
  specialised_outer<double> dbl_example2;
  specialised_outer<int> int_example2;

  outer& ptr2base1(dbl_example2);
  outer& ptr2base2(int_example2);

  //Check that results match the above
  ptr2base1.set_somevar(2);
  ptr2base1.Insides1.set_mA(0.4);
  ptr2base1.Insides2.set_mB(0.4);
  ptr2base2.set_somevar(2);
  ptr2base2.Insides1.set_mA(0.4);
  ptr2base2.Insides2.set_mB(0.4);

  std::cout << "From ptrs to base" << std::endl; 
  std::cout << ptr2base1.Insides1.get_f_of_mA(1.1) << std::endl;
  std::cout << ptr2base1.Insides2.get_f_of_mB(1.1, 1.1) << std::endl;
  std::cout << ptr2base2.Insides1.get_f_of_mA(1.1) << std::endl;
  std::cout << ptr2base2.Insides2.get_f_of_mB(1.1, 1.1) << std::endl;

}
