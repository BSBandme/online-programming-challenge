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


#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

ll dp[20][3][2];
struct LuckySum{
  long long construct(string note){
    reverse(note.begin(), note.end());
    int l = note.size();
    memset(dp, 0x3f, sizeof(dp));
    dp[0][2][0] = 0;
    for(ll i = 0, mi = 1; i < l; i++, mi *= 10) {
      int num = note[i] - '0';
      if(note[i] == '?') num = -1;
      for(int j = 1; j < 3; j++)
        for(int k = 0; k < 2; k++)
          if(dp[i][j][k] < 1ll << 60) {
            if(i != 0) {
              if(num == -1 || k + 4 == num)
                dp[i + 1][1][0] = min(dp[i][j][k] + (k + 4) * mi, dp[i + 1][1][0]);
              if(num == -1 || k + 7 == num)
                dp[i + 1][1][0] = min(dp[i][j][k] + (k + 7) * mi, dp[i + 1][1][0]);
            }
            if(j > 1) {
              if(num == -1 || k + 8 == num)
                dp[i + 1][2][0] = min(dp[i][j][k] + (k + 8) * mi, dp[i + 1][2][0]);
              if(num == -1 || k + 1 == num)
                dp[i + 1][2][1] = min(dp[i][j][k] + (k + 1) * mi, dp[i + 1][2][1]);
              if(num == -1 || k + 4 == num)
                dp[i + 1][2][1] = min(dp[i][j][k] + (k + 4) * mi, dp[i + 1][2][1]);
            }
          }
    }
    ll ans = min(dp[l][1][0], dp[l][2][0]);
    if(note[l - 1] == '?' || note[l - 1] == '1') {
      ll mi = 1;
      for(int i = 0; i < l - 1; i++)
        mi *= 10;
      ans = min(ans, dp[l - 1][2][1] + mi);
    }
    if(ans < 1ll << 60)
      return ans;

    return -1;
  }
  


}; 
