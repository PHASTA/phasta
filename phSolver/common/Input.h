#ifndef H_Input
#define H_Input
#include <vector>
#include <string>
#include <map>

#include "ValType.h"

using namespace std;

class Input {
public:
  Input(const string &, const string &default_fname = "");
  Input(const char*, const char* = "");
  ~Input();

  // return the entire input map
  map<string,string> InputMap() const;

  // returns the desired string
  //  const string &GetValue(const string &) const;
  ValType GetValue(const string &) const;

  // echo the entire input map
  void EchoInputMap(const ostream &ofile);

private:

  void trim_string(string *str);

  void get_input_lines(vector<string> *, ifstream& );
  void build_map(map<string,string> *, vector<string> *);

  map<string,string> *input_map;
  map<string,string> *default_map;

  vector<string> *input_text;
  vector<string> *default_text;

};




#endif
