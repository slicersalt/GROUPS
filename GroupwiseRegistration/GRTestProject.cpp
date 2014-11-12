#include <iostream>
//includes XML converted to C/C++ code
#include "CppTestProjectCLP.h"

using namespace std;

int main(int argc , char* argv[])
{
  //Should be the first line of your program to parse your command line
  PARSE_ARGS;
  //////////////////
  cout<<"Sphere: " << sphere << std::endl ;
  cout<<"tmpDepth: " << tmpDepth << std::endl ;
  cout<<"subjDepth: " << subjDepth << std::endl ;
  cout<<"coeff: " << coeff << std::endl ;
  cout<<"Correspondence: " << correspondence << std::endl ;
  cout<<"numSubj: " << numSubj << std::endl ;
  cout<<"Degree: " << degree << std::endl ;
  cout<<"Logfile: " << logfile << std::endl ;

/*
	-sphere = string1
	-tmpDepth = string2
	-subjDepth = list1 of strings
	-coeff = list2 of strings
	-correspondence = = list3 of strings (vector)
	-numSubj = integer (42)
	-degree = integer (9)
	-logfile = string3
	//GroupwiseRegistration(sphere, tmpDepth, subjDepth, coeff, correspondence, 42, 9, "log.txt");
*/
  return 0 ;
}
