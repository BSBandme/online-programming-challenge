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


#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

set <string> rs;

struct Bracket107{
  int yetanother(string s){
    rs.clear();
    int l = s.length();
    for(int i = 0; i < l; i++) {
      for(int j = 0; j < l; j++) if(s[j] != s[i]) {
        string str = s;
        char p = s[j];
        if(j < i) {
          for(int k = j; k < i; k++)
            str[k] = str[k + 1];
          str[i] = p;
        } else {
          for(int k = j; k > i; k--)
            str[k] = str[k - 1];
          str[i] = p;
        }
        bool flag = 1;
        for(int k = 0, cnt = 0; k < l; k++) {
          if(str[k] == '(') cnt++;
          else cnt--;
          if(cnt < 0) flag = 0;
        }
        if(flag)
          rs.insert(str);
      }
    }
    return rs.size();
  }
  


};

