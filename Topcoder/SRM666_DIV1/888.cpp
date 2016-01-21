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
	if(a >= mod) a -= mod;
}

#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

const int mod = MOD;

struct matrix{
	const static int maxn = 12;
	int row, col;
	ll mat[maxn][maxn];

	matrix(int r = 0, int c = 0){
		row = r; col = c;
		for(int i = 0; i < row; i++) for(int j = 0; j < col; j++) mat[i][j] = 0;
	}
	const ll * operator[] (int i) const {
		return mat[i];
	}
	ll * operator[] (int i) {
		return mat[i];
	}
	bool danweiju(){
		if(row != col) return 0;
		for(int i = 0; i < row; i++) for(int j = 0; j < col; j++) mat[i][j] = bool (i == j);
		return 1;
	}
	matrix operator * (const matrix& b) const{
		int i, j, k;
		matrix c(row, b.col);
		memset(c.mat, 0, sizeof(c.mat));
		for (i = 0; i < c.row; i++) for (k = 0; k < col; k++)
			if(mat[i][k])
				for (j = 0; j < c.col; j++){
					c.mat[i][j] += mat[i][k] * b.mat[k][j];
					c.mat[i][j] %= mod;
				}
		return c;
	}
	matrix operator + (const matrix& b) const{
		matrix c(max1(row, b.row), max1(col, b.col));
		for(int i = 0; i < c.row; i++) for(int j = 0; j < c.col; j++){
			ll a = 0; if(i < row && j < col) a = mat[i][j];
			ll b1 = 0; if(i < b.row && j < b.col) b1 = b.mat[i][j];
			c.mat[i][j] = a + b1;
			c.mat[i][j] %= mod;
		}
		return c;
	}
	inline void operator = (const matrix & b){
		memcpy(mat, b.mat, sizeof(mat));
		col = b.col;  row = b.row;
	}
	matrix pow(long long n){
		matrix ans(row, col), temp = *this;

		ans.danweiju();
		while(n){
			if(n & 1) ans = ans * temp;
			temp = temp * temp;
			n >>= 1;
		}
		return ans;
	}
};

int n, k, dn;
matrix  mat;

ll calc(int l, int r, int n, int cate) {
	if(l + r >= n) {
		if(l == n && r == n)
			if(n % 2 == cate)
				return 1;
		return 0;
	}
	int rn = n - l - r - 2;
	if(rn <= 0) {
		if((l + r) % 2 == cate)
			return 1;
		else
			return 0;
	}
	matrix  rmat = mat.pow(rn);
	ll ans = 0;
	cate ^= (l + r) % 2;
	for(int i = 0; i <= k; i++)
		add(ans, rmat[0][i + cate * (k + 1)]);
	return ans;
}
const int maxn = 55;
struct edge {
	int to, nxt;
} e[maxn * 2];
int head[maxn], le;
pii orig[maxn];
int fa[maxn];
ll dp[maxn][13][13][2];

void addedge(int a, int b) {
	e[le].to = b;
	e[le].nxt = head[a];
	head[a] = le++;
}

void justdoit(int no) {
	vector <pair <pii, int> > arr;
	for(int i = head[no]; i != -1; i = e[i].nxt)
		arr.push_back(mpr(orig[e[i].to], e[i].to));
	sort(arr.begin(), arr.end());
	int pn = arr.size();
	for(int i = 0; i < pn; i++)
		justdoit(arr[i].second);

	if(arr.size()) {
		for(int i = 0; i <= k; i++) for(int j = 0; j <= k; j++) {
			dp[no][i][j][0] = calc(i, j, arr[0].first.first - orig[no].first, 0);
			dp[no][i][j][1] = calc(i, j, arr[0].first.first - orig[no].first, 1);
		}
		int rl = arr[0].first.first - orig[no].first;
		for(int i = 0; i < (int)arr.size(); i++) {
			ll rdp[13][13][2];
			memset(rdp, 0, sizeof(rdp));
			int rno = arr[i].second;
			for(int ii = 0; ii <= k; ii++) for(int jj = 0; jj <= k; jj++) for(int cate = 0; cate < 2; cate++) if(dp[no][ii][jj][cate])
				for(int ii1 = 0; ii1 <= k - jj; ii1++) for(int jj1 = 0; jj1 <= k; jj1++) for(int rcate = 0; rcate < 1; rcate++) if(dp[rno][ii1][jj1][rcate]){
					if(jj1 == orig[rno].second - orig[rno].first + 1) {
						if(jj == rl)
							add(rdp[ii + ii1][jj + jj1][rcate ^ cate], dp[rno][ii1][jj1][rcate] * dp[no][ii][jj][cate] % mod);
						else
							add(rdp[ii][jj + jj1][rcate ^ cate], dp[rno][ii1][jj1][rcate] * dp[no][ii][jj][cate] % mod);
					} else {
						if(jj == rl)
							add(rdp[ii + ii1][jj1][rcate ^ cate], dp[rno][ii1][jj1][rcate] * dp[no][ii][jj][cate] % mod);
						else
							add(rdp[ii][jj1][rcate ^ cate], dp[rno][ii1][jj1][rcate] * dp[no][ii][jj][cate] % mod);
					}
				}
			memcpy(dp[no], rdp, sizeof(rdp));
			memset(rdp, 0, sizeof(rdp));
			rl += orig[rno].second - orig[rno].first + 1;

			int tl;
			if(i != (int)arr.size() - 1)
				tl = arr[i + 1].first.first - arr[i].first.second - 1;
			else tl = orig[no].second - arr[i].first.second;
			if(tl == 0) continue;
			for(int ii1 = 0; ii1 <= k; ii1++) for(int jj1 = 0; jj1 <= k; jj1++) for(int rcate = 0; rcate < 2; rcate++){
				int tdp = calc(ii1, jj1, tl, rcate);
				if(!tdp) continue;
				for(int ii = 0; ii <= k; ii++) for(int jj = 0; jj <= k - ii1; jj++) for(int cate = 0; cate < 2; cate++) if(dp[no][ii][jj][cate]) {

					if(jj1 == tl) {
						if(jj == rl)
							add(rdp[ii + jj1][jj + jj1][rcate ^ cate], tdp * dp[no][ii][jj][cate] % mod);
						else
							add(rdp[ii][jj + jj1][rcate ^ cate], tdp * dp[no][ii][jj][cate] % mod);
					} else
						if(jj == rl)
							add(rdp[ii + ii1][jj1][rcate ^ cate], tdp * dp[no][ii][jj][cate] % mod);
						else
							add(rdp[ii][jj1][rcate ^ cate], tdp * dp[no][ii][jj][cate] % mod);
				}
			}
			memcpy(dp[no], rdp, sizeof(rdp));
			rl += tl;
		}
	} else
		for(int i = 0; i <= k; i++) for(int j = 0; j <= k; j++) for(int cate = 0; cate < 2; cate++)
			dp[no][i][j][cate] = calc(i, j, orig[no].second - orig[no].first + 1, cate);
}

struct CountBinarySequences{
	int countSequences(int rn, int rk, vector <int> L, vector <int> R){
		memset(mat.mat, 0, sizeof(mat.mat));
		n = rn, k = rk;
		mat.col = mat.row = (k + 1) * 2;
		for(int i = 0; i < k; i++) {
			mat[i][0] = 1;
			mat[i][i + 1 + k + 1] = 1;
			mat[i + k + 1][k + 1] = 1;
			mat[i + k + 1][i + 1] = 1;
		}
		mat[k][0] = 1;
		mat[k + k + 1][k + 1] = 1;
		memset(head, -1, sizeof(head));
		dn = L.size();
		for(int i = 0; i < dn; i++)
			orig[i] = mpr(L[i] - 1, R[i] - 1);
		orig[dn] = mpr(0, n - 1);
		for(int i = 0; i < dn; i++) {
			fa[i] = dn;
			for(int j = 0; j < dn; j++) if(j != i) {
				if(orig[j].first <= orig[i].first && orig[j].second >= orig[i].second) {
					if(fa[i] != dn && orig[j].second - orig[j].first > orig[fa[i]].second - orig[fa[i]].first)
						continue;
					fa[i] = j;
				}
			}
		}
		for(int i = 0; i < dn; i++)
			addedge(fa[i], i);

		justdoit(dn);
		ll ans = 0;
		for(int i = 0; i <= k; i++) for(int j = 0; j <= k; j++)
			for(int tt = 0; tt < 2; tt++)
				add(ans, dp[dn][i][j][tt]);

		return int(ans) ;
	}



};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

