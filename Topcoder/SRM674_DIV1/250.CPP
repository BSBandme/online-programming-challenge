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


struct VampireTree{
  int maxDistance(vector <int> num){
    int n = num.size();
    int cnt = 0;
    for(int i = 0; i < n; i++)
      cnt += num[i] - 1;
    cnt++;
    if(cnt != n - 1)
      return -1;

    int ans = 0;
    for(int i = 0; i < n; i++)
      if(num[i] != 1)
        ans++;

    return ans + 1;
  }
  


}; 
