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

long long pow(long long n, long long m, long long mod = 0){
	if(m < 0) return 0;
	long long ans = 1;
	long long k = n;
	while(m){
		if(m & 1) {
			ans *= k;
			if(mod) ans %= mod;
		}
		k *= k;
		if(mod) k %= mod;
		m >>= 1;
	}
	return ans;
}

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
const int maxn = 100;
ll mi[maxn];
ll xi[maxn][maxn], c[maxn][maxn], bit[maxn];
ll ans[maxn][maxn];

struct PowerPartition{
	int count(int m, long long x){
		xi[1][2] = xi[1][1] = (mod + 1) / 2;
		for(int i = 2; i < maxn - 2; i++) {
			for(int j = 2; j <= i + 1; j++)
				xi[i][j] = xi[i - 1][j - 1] * i % mod * pow(j, mod - 2, mod) % mod;
			xi[i][1] = 1;
			for(int j = 2; j <= i + 1; j++) {
				xi[i][1] += mod - xi[i][j];
				xi[i][1] %= mod;
			}
		}
		for(int i = 0; i < maxn; i++) {
			c[i][0] = 1;
			for(int j = 1; j <= i; j++) {
				c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
				c[i][j] %= mod;
			}
		}
		memset(bit, 0, sizeof(bit));
		for(ll i = 0, rx = x; rx; i++, rx /= m)
			bit[i] = rx % m;

		memset(ans, 0, sizeof(ans));
		ans[0][0] = 1;
		for(ll i = 1; i < 70; i++) {
			ll temp[maxn];
			memset(temp, 0, sizeof(temp));
			for(int j = 0; j <= i; j++) {
				for(int k = 0; k <= j; k++)
					add(temp[k], c[j][k] * pow(m, k, mod) % mod * pow(bit[i - 1], j - k, mod) % mod * ans[i - 1][j] % mod);
			}
			ans[i][1] = temp[0], ans[i][0] = temp[0];
			for(int j = 1; j <= i; j++) {
				for(int k = 0; k <= j + 1; k++)
					add(ans[i][k], xi[j][k] * temp[j] % mod);
			}
		}

		return ans[69][0] ;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

