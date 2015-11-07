#ifndef H_Input
#define H_Input
#include <vector>
#include <string>
#include <map>

#include "ValType.h"

using namespace std;

namespace phSolver{
class Input {
public:
  Input(const string &, const string &default_fname = "");
  ~Input();

  // return the entire input map
  map<string,string> InputMap() const;

  // returns the desired string
  //  const string &GetValue(const string &) const;
  ValType GetValue(const string &) const;

  // echo the entire input map
  void EchoInputMap(const ostream &ofile);

  const char* GetUserFileName();
  const char* GetDefaultFileName();
private:

  void trim_string(string *str);

  void get_input_lines(vector<string> *, ifstream& );
  void build_map(map<string,string> *, vector<string> *);

  map<string,string> *input_map;
  map<string,string> *default_map;

  vector<string> *input_text;
  vector<string> *default_text;

  string userConfFileName;
  string defaultConfFileName;

};
} //end phSolver namespace

#endif
