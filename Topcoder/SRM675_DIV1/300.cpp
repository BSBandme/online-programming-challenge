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

typedef vector <int> vi;


#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................


struct TreeAndPathLength3{

  vector <int> construct(int s){
    int rn = 2;
    for(; rn * (rn - 1) <= s; rn++);
    rn--;
    vi ans;
    for(int i = 0; i < rn; i++) {
      ans.push_back(0);
      ans.push_back(i * 2 + 1);
      ans.push_back(i * 2 + 1);
      ans.push_back(i * 2 + 2);
    }
    for(int j = rn * (rn - 1), i = rn * 2; j < s; j++, i++) {
      ans.push_back(i);
      ans.push_back(i + 1);
    }

    return ans ;
  }
  


}; 
