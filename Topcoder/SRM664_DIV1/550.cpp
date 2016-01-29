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
const int maxn = 1001000;
struct edge {
	int to, nxt;
} e[maxn * 2];
int head[maxn], le;
ll fac[maxn], inv[maxn], invfac[maxn];
ll expc[maxn];
int n;
ll ans;

inline void addedge(int a, int b) {
	e[le].to = b;
	e[le].nxt = head[a];
	head[a] = le++;
}
int dfs(int no, int fa) {
	expc[no] = 0;
	for(int i = head[no]; i != -1; i = e[i].nxt) if(e[i].to != fa) {
		add(expc[no], dfs(e[i].to, no));
	}
	expc[no] = (expc[no] + 1) * inv[no + 1] % mod;
	return expc[no];
}
void dfs1(int no, int fa = -1, ll already = 0) {
	add(ans, (2 * already + 1) * inv[no + 1] % mod);
	ll texp = already;
	for(int i = head[no]; i != -1; i = e[i].nxt) if(e[i].to != fa){
		dfs1(e[i].to, no, (texp + 1) * inv[no + 1] % mod);
		add(texp, expc[e[i].to]);
	}
}
struct BearAttacks{
	int expectedValue(int N, int R0, int A, int B, int M, int LOW, int HIGH){
		memset(head, -1, sizeof(head));
		le = 0;
		n = N;

		ll R = R0;
		ll BILLION = 1000000000;
		for(int i = 1; i < n; i++) {
		    R = (A * R + B) % M;
		    ll MIN = (1LL * (i-1) * LOW)  / BILLION;
		    ll MAX = (1LL * (i-1) * HIGH) / BILLION;
		    ll no = MIN + (R % (MAX-MIN+1));
		    addedge(i, no);
		    addedge(no, i);
		}

		invfac[0] = fac[0] = 1;
		for(int i = 1; i < maxn; i++) {
			fac[i] = fac[i - 1] * i % mod;
			inv[i] = pow(i, mod - 2, mod);
			invfac[i] = invfac[i - 1] * inv[i] % mod;
		}

		ans = 0;
		dfs(0, -1);
		dfs1(0, -1);
		ans *= fac[n];
		ans %= mod;


		return ans;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor


// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

