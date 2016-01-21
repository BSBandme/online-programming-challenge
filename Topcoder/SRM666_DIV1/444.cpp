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

const int mod = MOD;
const int maxn = 4010;
ll dp[maxn][2][2], c[maxn][maxn];
int rn;

struct SumOverPermutations{
	ll justdoit(int n, int fl, int fr) {
		if(n == 0) return 1;
		if(n == 1)
			return rn - fl - fr;
		if(dp[n][fl][fr] != -1) return dp[n][fl][fr];
		dp[n][fl][fr] = 0;

		for(int i = 1; i < n - 1; i++) {
			add(dp[n][fl][fr], rn * justdoit(i, fl, 1) % mod * justdoit(n - 1 - i, 1, fr) % mod * c[n - 1][i] % mod);
		}
		add(dp[n][fl][fr], (rn - fl) * justdoit(n - 1, 1, fr) % mod);
		add(dp[n][fl][fr], (rn - fr) * justdoit(n - 1, fl, 1) % mod);
		return dp[n][fl][fr];
	}
	int findSum(int n){
		rn = n;
		for(int i = 0; i < maxn; i++) {
			c[i][0] = 1;
			for(int j = 1; j <= i; j++) {
				c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
				c[i][j] %= mod;
			}
		}
		memset(dp, -1, sizeof(dp));
		return justdoit(n, 0, 0);
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

