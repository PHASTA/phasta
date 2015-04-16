#ifndef H_ValType
#define H_ValType

#include <iostream>
#include <string>

using namespace std;

class ValType {
public:
  ValType(const string & s)
    : str(s)
    { used = false; }

  ~ValType()
    {
      if (!used) 
	cerr << "error: ambiguous return type" << endl;
    }

  // conversion operators
  operator double();
  operator double*();
  operator vector<double>();
  operator vector<int>();
  operator int();
  operator string();

private:
  bool used;
  string str;

  int get_int(string str);
  double get_double(string str);
  double *get_double_array(string str);
  vector<double> get_vector(string str);
  vector<int> get_ivector(string str);
  string get_string(string str);
  
};


#endif

