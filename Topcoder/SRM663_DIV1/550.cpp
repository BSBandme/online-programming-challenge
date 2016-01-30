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

typedef long long type;
struct comp{
    double x, y;
    comp(double _x=0, double _y=0) : x(_x), y(_y) {}
};
namespace FFT{
    const int N = 1<<12;
    const double pi2 = 3.1415926535897932 * 2;
    comp a[N], b[N], tmp[N];
    int n, bn;
    type res[N];
    inline comp W(int n, bool inv) {
        double ang = inv ? -pi2 / n : pi2 / n;
        return comp(cos(ang), sin(ang));
    }
    int bitrev(int x) {
        int ans = 0;
        for (int i=1; i<=bn; ++i)
            ans <<= 1, ans |= x & 1, x >>= 1;
        return ans;
    }
    void dft(comp *a,bool inv) {
        int step, to; comp w, wi, A, B;
        for (int i=0; i<n; ++i) {
            to = bitrev(i);
            if (to > i) std::swap(a[to], a[i]);
        }
        for (int i=1; i<=bn; ++i) {
            wi = W(1<<i, inv); w = comp(1, 0);
            step = 1 << (i-1);
            for (int k=0; k<step; ++k) {
                for (int j=0; j<n; j+=1<<i) {
                    int t = j | k, d = j|k|step;
                    A = a[t];
                    B.x  = w.x * a[d].x - w.y * a[d].y;
                    B.y  = w.x * a[d].y + w.y * a[d].x;
                    a[t].x = A.x + B.x, a[t].y = A.y + B.y;
                    a[d].x = A.x - B.x, a[d].y = A.y - B.y;
                }
                comp tmp;
                tmp.x = w.x * wi.x - w.y * wi.y;
                tmp.y = w.x * wi.y + w.y * wi.x;
                w = tmp;
            }
        }
    }
    int mul(int n1, ll *x1, int n2, ll *x2) {
        n = std::max(n1, n2);
        for (bn = 0; (1<<bn) < n; ++bn); ++bn;
        n = 1 << bn;
        for (int i=0; i<n; ++i) a[i] = b[i] = comp(0, 0);
        for (int i=0; i<n1; ++i) a[i] = comp(x1[i], 0);
        for (int i=0; i<n2; ++i) b[i] = comp(x2[i], 0);
        dft(a, false); dft(b, false);
        for (int i=0; i<n; ++i) {
            tmp[i].x = a[i].x * b[i].x - a[i].y * b[i].y;
            tmp[i].y = a[i].x * b[i].y + a[i].y * b[i].x;
        }
        dft(tmp, true);
        for (int i=0; i<n; ++i) res[i] = (type)(tmp[i].x/n + 0.3);
        for (--n; n && !res[n]; --n);
        return n+1;
    }
}

const int mod = MOD;
const int maxn = 2010;
const int rmaxn = 1000100;
int minc[maxn];
int arr[maxn];
int d, rn;
ll dp[2][maxn][maxn];

ll inv[rmaxn], fac[rmaxn], invfac[rmaxn], tfac[rmaxn], tinvfac[rmaxn];

int c(int a, int b) {
	a %= mod;
	if(a < b) return 0;
	if(a < rmaxn)
		return fac[a] * invfac[b] % mod * invfac[a - b] % mod;
	return tfac[mod - a + b - 1] * tinvfac[mod - a - 1] % mod * invfac[b] % mod;
}

struct ChangingChange{
	vector <int> countWays(vector <int> ways, vector <int> vr, vector <int> nr){
		memset(minc, 0, sizeof(minc));
		d = ways.size() - 1;
		rn = vr.size();
		tfac[0] = tinvfac[0] = fac[0] = invfac[0] = 1;
		for(int i = 1; i < rmaxn; i++) {
			inv[i] = pow(i, mod - 2, mod);
			invfac[i] = invfac[i - 1] * inv[i] % mod;
			fac[i] = fac[i - 1] * i % mod;
			tfac[i] = tfac[i - 1] * (mod - i) % mod;
			tinvfac[i] = pow(tfac[i], mod - 2, mod);
		}
		memset(minc, 0, sizeof(minc));
		for(int i = 0; i < rn; i++) {
			minc[vr[i]] = max(minc[vr[i]], nr[i]);
		}
		memset(dp, 0, sizeof(dp));
		dp[0][0][0] = 1;
		for(int i = 1; i <= d; i++) {
			arr[i] = ways[i] - dp[0][i - 1][i];
			if(arr[i] < minc[i])
				arr[i] += mod;
			memcpy(dp[0][i], dp[0][i - 1], sizeof(dp[0][i]));
			for(int j = i; j <= d; j++) {
				ll c = 1;
				for(int k = 1; k * i <= j; k++) {
					c = c * (arr[i] - k + 1) % mod * inv[k] % mod;
					add(dp[0][i][j], dp[0][i - 1][j - k * i] * c % mod);
				}
			}
		}
		dp[1][d + 1][0] = 1;
		for(int i = d; i >= 1; i--) {
			FFT::mul(d + 1, dp[0][i - 1], d + 1, dp[1][i + 1]);
			for(int j = 0; j <= d; j++)
				dp[0][i][j] = FFT::res[j] % mod;
			memcpy(dp[1][i], dp[1][i + 1], sizeof(dp[0][i]));
			for(int j = i; j <= d; j++) {
				ll c = 1;
				for(int k = 1; k * i <= j; k++) {
					c = c * (arr[i] - k + 1) % mod * inv[k] % mod;
					add(dp[1][i][j], dp[1][i + 1][j - k * i] * c % mod);
				}
			}
		}
		vi rans;
		for(int i = 0; i < rn; i++) {
			int num = vr[i], all = arr[num] - nr[i];
			ll ans = dp[0][num][d];
			for(int j = 1; j * num <= d; j++) {
				ans += dp[0][num][d - j * num] * c(all, j) % mod;
				ans %= mod;
			}
			rans.push_back(ans);
		}

		return rans;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

