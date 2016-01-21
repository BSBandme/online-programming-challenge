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

#define mpr make_pair
typedef unsigned int ui;
typedef unsigned long long ull;
typedef long long ll;

template <class T> inline T abs1(T a) {return a < 0 ? -a : a;}

const int maxn = 210;
ll dp[maxn][maxn];
ll fac[maxn];
ll c[maxn][maxn];
ll pw[maxn][maxn * maxn];
const int mod = 1000000007;
int rl;

ll getdp(int n, int l) {
  if(n == 0) return 1;
  if(l == 0) return 0;
  if(n == 1) return l;
  if(dp[n][l] != -1)
    return dp[n][l];
  dp[n][l] = 0;
  for(int i = 1; i <= l; i++)
    for(int j = 0; j < n; j++) {
      dp[n][l] += getdp(j, i - 1) * getdp(n - j - 1, i) % mod * pw[rl - i + 1][(n - j) * (j + 1) - 1];
      dp[n][l] %= mod;
    }
  return dp[n][l];
}

struct LineMST{
  int count(int n, int l){
    for(int i = 0; i < maxn; i++) {
      c[i][0] = 1;
      for(int j = 1; j <= i; j++) {
        c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
        c[i][j] %= mod;
      }
    }
    fac[0] = 1;
    for(int j = 1; j < maxn; j++) {
      fac[j] = fac[j - 1] * j % mod;
    }
    for(int i = 0; i < maxn; i++) {
      pw[i][0] = 1;
      for(int j = 1; j < maxn * maxn; j++)
        pw[i][j] = pw[i][j - 1] * i % mod;
    }
    rl = l;
    memset(dp, 0, sizeof(dp));
    for(int i = 0; i <= l; i++)
      dp[0][i] = 1;
    for(int i = 1; i <= l; i++)
      dp[1][i] = i;
    for(int ii = 2; ii < n; ii++) {
      for(int jj = 1; jj <= l; jj++) {
        dp[ii][jj] = dp[ii][jj - 1];
          for(int j = 0; j < ii; j++) {
            dp[ii][jj] += dp[j][jj - 1] * dp[ii - j - 1][jj] % mod * pw[l - jj + 1][(ii - j) * (j + 1) - 1];
            dp[ii][jj] %= mod;
          }
      }
    }

    return dp[n - 1][l] * fac[n] % mod * ((mod + 1) / 2) % mod;
  }
  


}; 
