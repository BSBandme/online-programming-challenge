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


#define  MOD 1000000007
template <class t1, class t2>
inline void add(t1 &a, t2 b, int mod = -1) {
  if(mod == -1) mod = MOD;
  a += b;
  while(a >= mod) a -= mod;
  while(a < 0) a += mod;
}

#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

const int maxn = 55;
int dp[maxn][maxn * 2];

struct ExactTree{
  int getTree(int n, int m, int r){
    memset(dp, 0x3f, sizeof(dp));
    dp[1][0] = 0;
    for(int i = 2; i <= n; i++)  {
      for(int j = 1; j < i; j++) for(int m1 = 0; m1 < m; m1++) if(dp[j][m1] < MOD) {
        for(int m2 = 0; m2 < m; m2++) if(dp[i - j][m2] < MOD)
          dp[i][(m1 + m2 + j * (n - j)) % m] = min(dp[i][(m1 + m2 + j * (n - j)) % m], dp[j][m1] + dp[i - j][m2] + j * (n - j));
      }
    }
    if(dp[n][r] > MOD)
      return -1;
    return dp[n][r] ;
  }
  


}; 
