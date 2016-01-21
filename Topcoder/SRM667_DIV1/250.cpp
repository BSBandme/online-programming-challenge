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
typedef pair <int, int> pii;
typedef pair <ll, ll> pll;
typedef pair <double, double> pdd;
typedef vector <int> vi;
typedef vector <ll> vll;
typedef vector <double> vd;
typedef vector <string> vs;
typedef map <string, int> mpsi;
typedef map <double, int> mpdi;
typedef map <int, int> mpii;


template <class T> inline T lowb(T num) {return num & (-num); }
inline int bitnum(ui nValue) { return __builtin_popcount(nValue); }
inline int bitnum(int nValue) { return __builtin_popcount(nValue); }
inline int bitnum(ull nValue) { return __builtin_popcount(nValue) + __builtin_popcount(nValue >> 32); }
inline int bitnum(ll nValue) { return __builtin_popcount(nValue) + __builtin_popcount(nValue >> 32); }
inline int bitmaxl(ui a) { if(a == 0) return 0; return 32 - __builtin_clz(a); }
inline int bitmaxl(int a) { if(a == 0) return 0; return 32 - __builtin_clz(a); }
inline int bitmaxl(ull a) { int temp = a >> 32; if(temp) return 32 - __builtin_clz(temp) + 32; return bitmaxl(int(a)); }
inline int bitmaxl(ll a) { int temp = a >> 32; if(temp) return 32 - __builtin_clz(temp) + 32; return bitmaxl(int(a)); }

#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

const int maxn = 55;
int arr[maxn];
int dp[1 << 21];
int n, l;

struct OrderOfOperations{
  int minTime(vector <string> s){
    n = s.size();
    l = s[0].size();
    int all = 0;
    for(int i = 0; i < n; i++) {
      arr[i] = 0;
      for(int j = 0; j < l; j++) {
        arr[i] = arr[i] * 2 + int(s[i][j] - '0');
      }
    }
    for(int i = 0; i < n; i++)
      all |= arr[i];

    memset(dp, 0x3f, sizeof(dp));
    dp[0] = 0;
    for(int j = 0; j < 1 << l; j++) if((j & all) == j){
      for(int i = 0; i < n; i++) if((j & arr[i]) != arr[i]) {
        int temp = arr[i] ^ (arr[i] & j);
        int cnt = bitnum(temp);
        dp[j | arr[i]] = min(dp[j | arr[i]], dp[j] + cnt * cnt);
      }
    }

    return dp[all] ;
  }
  


};

