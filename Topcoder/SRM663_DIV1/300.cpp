/*
 * temp.cpp
 *
 *  Created on: 2012-7-18
 *      Author: BSBandme
 */
//#pragma comment(linker, "/STACK:1024000000,1024000000")
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <cassert>
#include <list>
#include <iomanip>
#include <math.h>
#include <deque>
#include <utility>
#include <map>
#include <set>
#include <bitset>
#include <numeric>
#include <climits>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <sstream>
#include <tr1/unordered_set>
#include <tr1/unordered_map>

using namespace std;
using namespace tr1;


string rini;

struct ABBADiv1{
  int can(string str) {
    if(str.length() == rini.length())
      return rini == str;
    if(str[0] == 'B') {
      string tstr = str.substr(1, str.size() - 1);
      reverse(tstr.begin(), tstr.end());
      if(can(tstr)) return 1;
    }
    if(str[str.size() - 1] == 'A') {
      string tstr = str.substr(0, str.size() - 1);
      return can(tstr);
    }
    return 0;

  }
  string canObtain(string initial, string target){
    rini = initial;
    if(can(target)) return "Possible";
    return "Impossible";
  }
  


};

