#include <string>
#include <vector>
//MR CHANGE
#include <cstring>
#include <cmath>
#include <cstdlib>
//MR CHANGE END

#include "ValType.h"

ValType::operator int()
{
  used = true;
  return get_int(str);
}

ValType::operator vector<double>()
{
  used = true;
  return get_vector(str);
}

ValType::operator vector<int>()
{
  used = true;
  return get_ivector(str);
}

ValType::operator double()
{
  used = true;
  return get_double(str);
}

ValType::operator double*()
{
  used = true;
  return get_double_array(str);
}

ValType::operator string()
{
  used = true;
  return get_string(str);
}

//
//  function implementations for type specific conversions
//

int ValType::get_int(string str)
{
  int i = atoi(str.c_str());
  return i;
}

double ValType::get_double(string str)
{
  double x = atof(str.c_str());
  return x;
}

double *ValType::get_double_array(string str)
{
  //istrstream ist(str.c_str(),str.length());
  //vector<double> vec;
  //double v;
  //while ( ist >> v ) {
  //  vec.push_back(v);
  //}
  vector<double> vec = get_vector(str);
  int n = vec.size();
  double *x = new double[n];
  for (int i=0; i < n; i++) {
    x[i] = vec[i];
  }
  return x;
}

vector<double> ValType::get_vector(string str)
{
  char* s = (char*) malloc(str.size()+1);
  strcpy(s,str.c_str());
  char *strTok = strtok(s," ");
  vector<double> vec;
  while(strTok) {
    vec.push_back(atof(strTok));
    strTok = strtok(NULL," ");
  } 
  free(s);
  return vec;
}

vector<int> ValType::get_ivector(string str)
{
  char* s = (char*) malloc(str.size()+1);
  strcpy(s,str.c_str());
  char *strTok = strtok(s," ");
  vector<int> vec;
  while(strTok) {
    vec.push_back(atoi(strTok));
    strTok = strtok(NULL," ");
  }
  free(s);
  return vec;
}

string ValType::get_string(string str)
{
  return str;
}
