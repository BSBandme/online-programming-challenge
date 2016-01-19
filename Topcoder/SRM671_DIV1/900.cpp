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

const double pi = acos(0.0) * 2.0;
const double eps = 1e-12;
const int step[8][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

template <class T> inline T abs1(T a) {return a < 0 ? -a : a;}

template <class T> inline T max1(T a, T b) { return a > b ? a : b; }
template <class T> inline T max1(T a, T b, T c) { return max1(max1(a, b), c); }
template <class T> inline T max1(T a, T b, T c, T d) { return max1(max1(a, b, c), d); }
template <class T> inline T max1(T a, T b, T c, T d, T e) { return max1(max1(a, b, c, d), e); }
template <class T> inline T min1(T a, T b) { return a < b ? a : b; }
template <class T> inline T min1(T a, T b, T c) { return min1(min1(a, b), c); }
template <class T> inline T min1(T a, T b, T c, T d) { return min1(min1(a, b, c), d); }
template <class T> inline T min1(T a, T b, T c, T d, T e) { return min1(min1(a, b, c, d), e); }

inline int jud(double a, double b){
	if(abs(a) < eps && abs(b) < eps) return 0;
	else if(abs1(a - b) / abs1(a) < eps) return 0;
	if(a < b) return -1;
	return 1;
}
template <typename t> inline int jud(t a, t b){
	if(a < b) return -1;
	if(a == b) return 0;
	return 1;
}

// f_lb == 1代表返回相同的一串的左边界，f_small == 1代表返回如果没有寻找的值返回小的数
template <typename it, typename t1>
inline int find(t1 val, it a, int na, bool f_small = 1, bool f_lb = 1){
	int be = 0, en = na - 1;
	if(*a <= *(a + na - 1)){
		if(f_lb == 0) while(be < en){
			int mid = (be + en + 1) / 2;
			if(jud(*(a + mid), val) != 1) be = mid;
			else en = mid - 1;
		}else while(be < en){
			int mid = (be + en) / 2;
			if(jud(*(a + mid), val) != -1) en = mid;
			else be = mid + 1;
		}
		if(f_small && jud(*(a + be), val) == 1) be--;
		if(!f_small && jud(*(a + be), val) == -1) be++;
	} else {
		if(f_lb) while(be < en){
			int mid = (be + en + 1) / 2;
			if(jud(*(a + mid), val) != -1) be = mid;
			else en = mid - 1;
		}else while(be < en){
			int mid = (be + en) / 2;
			if(jud(*(a + mid), val) != 1) en = mid;
			else be = mid + 1;
		}
		if(!f_small && jud(*(a + be), val) == -1) be--;
		if(f_small && jud(*(a + be), val) == 1) be++;
	}
	return be;
}

template <class T> inline T lowb(T num) {return num & (-num); }
inline int bitnum(ui nValue) { return __builtin_popcount(nValue); }
inline int bitnum(int nValue) { return __builtin_popcount(nValue); }
inline int bitnum(ull nValue) { return __builtin_popcount(nValue) + __builtin_popcount(nValue >> 32); }
inline int bitnum(ll nValue) { return __builtin_popcount(nValue) + __builtin_popcount(nValue >> 32); }
inline int bitmaxl(ui a) { if(a == 0) return 0; return 32 - __builtin_clz(a); }
inline int bitmaxl(int a) { if(a == 0) return 0; return 32 - __builtin_clz(a); }
inline int bitmaxl(ull a) { int temp = a >> 32; if(temp) return 32 - __builtin_clz(temp) + 32; return bitmaxl(int(a)); }
inline int bitmaxl(ll a) { int temp = a >> 32; if(temp) return 32 - __builtin_clz(temp) + 32; return bitmaxl(int(a)); }

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
template <class t>
void output1(t arr) {
	for(int i = 0; i < (int)arr.size(); i++)
		cerr << arr[i] << ' ';
	cerr << endl;
}
template <class t>
void output2(t arr) {
	for(int i = 0; i < (int)arr.size(); i++)
		output1(arr[i]);
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

const int maxn = 33;
ll dp[maxn][15][1 << 14][2];
int p[100][5];

struct BearDestroys{
	int sumUp(int w, int h, int mod){
		memset(p, 0, sizeof(p));
		for(int i = 0; i < w + h - 1; i++) {
			int stw = min(i, w - 1), sth = i - stw;
			int enh = min(i, h - 1), enw = i - enh;
			p[i][0] = stw, p[i][1] = sth, p[i][2] = enw, p[i][3] = enh;
			if(i) {
				if(p[i][0] > p[i - 1][0])
					p[i - 1][4] |= 1;
				if(p[i][3] > p[i - 1][3])
					p[i - 1][4] |= 2;
			}
		}
		memset(dp, 0, sizeof(dp));
		dp[0][0][0][0] = 1;
		dp[0][0][0][1] = 0;
		for(int i = 0; i < w + h - 2; i++) {
			int stw = p[i][0], sth = p[i][1], enw = p[i][2], enh = p[i][3], rf = p[i][4];
			if(rf & 1) {
				for(int j = (1 << (enh - sth + 1)) - 1; j >= 0; j--) {
					dp[stw][sth][j << 1][0] = dp[stw][sth][j][0];
					dp[stw][sth][j << 1][1] = dp[stw][sth][j][1];
					if(j & 1) {
						dp[stw][sth][j][0] = 0;
						dp[stw][sth][j][1] = 0;
					}
				}
			}
			for(int rw = stw, rh = sth, j = rf & 1; rh <= enh; rh++, rw--, j++) {
				int nxtj = rw - 1, nxti = rh + 1;
				if(rh == enh) {
					nxtj = p[i + 1][0];
					nxti = p[i + 1][1];
				}
				for(int mask = 0; mask < 1 << (enh - sth + 1 + (rf & 1)); mask++) {
					if(mask & (1 << j)) {
						add(dp[nxtj][nxti][mask ^ (1 << j)][0], dp[rw][rh][mask][0] * 2, mod);
						add(dp[nxtj][nxti][mask ^ (1 << j)][1], dp[rw][rh][mask][1] * 2, mod);
					} else {
						if(j && (mask & (1 << (j - 1))) == 0) {
							add(dp[nxtj][nxti][mask ^ (1 << (j - 1))][0], dp[rw][rh][mask][0], mod);
							add(dp[nxtj][nxti][mask ^ (1 << (j - 1))][1], dp[rw][rh][mask][1] + dp[rw][rh][mask][0], mod);
						} else {
							if(bool(rf & 2) || (rw != enw)) {
								add(dp[nxtj][nxti][mask | (1 << j)][0], dp[rw][rh][mask][0], mod);
								add(dp[nxtj][nxti][mask | (1 << j)][1], dp[rw][rh][mask][1] + dp[rw][rh][mask][0], mod);
							} else {
								add(dp[nxtj][nxti][mask][0], dp[rw][rh][mask][0], mod);
								add(dp[nxtj][nxti][mask][1], dp[rw][rh][mask][1], mod);
							}
						}
						if(bool(rf & 2) || (rw != enw)) {
							add(dp[nxtj][nxti][mask | (1 << j)][0], dp[rw][rh][mask][0], mod);
							add(dp[nxtj][nxti][mask | (1 << j)][1], dp[rw][rh][mask][0] + dp[rw][rh][mask][1], mod);
						} else {
							if(mask & (1 << (j - 1))) {
								add(dp[nxtj][nxti][mask][0], dp[rw][rh][mask][0], mod);
								add(dp[nxtj][nxti][mask][1], dp[rw][rh][mask][1], mod);
							} else {
								add(dp[nxtj][nxti][mask ^ (1 << (j - 1))][0], dp[rw][rh][mask][0], mod);
								add(dp[nxtj][nxti][mask ^ (1 << (j - 1))][1], dp[rw][rh][mask][1] + dp[rw][rh][mask][0], mod);
							}
						}
					}
				}
			}
		}

		return (dp[w - 1][h - 1][0][1] + dp[w - 1][h - 1][1][1]) * 2 % mod;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

