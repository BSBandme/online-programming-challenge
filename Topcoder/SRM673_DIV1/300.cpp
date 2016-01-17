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
int rw[maxn];
int rh[maxn];
int m;
const int mod = MOD;

struct BearCavalry{
  int countAssignments(vector <int> w, vector <int> h){
    int n = w.size();
    int rn = n - 1;
    for(int j = 1; j < n; j++)
      rw[j - 1] = w[j];
    sort(rw, rw + n - 1);
    int rans = 0;
    for(int i = 0; i < n; i++) {
      int lh = 0;
      for(int j = 0; j < n; j++) if(j != i)
        rh[lh++] = h[j];
      sort(rh, rh + n - 1);
      int maxv = w[0] * h[i];
      ll ans = 1;
      int rp = 0;
      for(int j = n - 2; j >= 0; j--) {
        while(rh[rp] * rw[j] < maxv && rp < n - 1)
          rp++;
        ans = ans * (rp - (n - j) + 2) % mod;
      }
      add(rans, ans);
    }
    return rans;
  }
  


};


