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

const int maxv = 20000;
ll dp[maxv + 20000];
ll dp1[maxv + 20000];
struct pac{
	ll cnt, w, v;
	bool operator < (const pac &a) const {
		return a.v * w > v * a.w;
	}
	pac(int a = 0, int b = 0, int c = 0) :
		cnt(a), w(b), v(c) {}
} arr[100];

multiset <ll> ms;

void getdp(ll dp[maxv], int cate, vector <pac> &rarr) {
	int rn = rarr.size();
	for(int i = 0; i < rn; i++) {
		ms.clear();
		ll rv = rarr[i].v, rcnt = rarr[i].cnt, rw = rarr[i].w;
		for(int mod = 0; mod < rw; mod++) {
			ll st = maxv / rw * rw + mod;
			if(st > maxv) st -= rw;
			ll stc = st / rw;
			for(int j = 0; j <= rcnt && j * rw <= st; j++)
				ms.insert(dp[st - j * rw] + rv * j);
			for(; st >= 0; st -= rw) {
				ll rdp = dp[st] + (rv * (stc - st / rw));
				if(cate == 1)
					dp[st] = max(dp[st], *ms.rbegin() - rv * (stc - st / rw));
				else
					dp[st] = min(dp[st], *ms.begin() - rv * (stc - st / rw));
				auto it = ms.find(rdp);
				assert(it != ms.end());
				ms.erase(it);
//				cerr << ms.size() << endl;
				if(st - (rcnt + 1) * rw >= 0)
					ms.insert(dp[st - (rcnt + 1) * rw] + rv * (stc - (st - (rcnt + 1) * rw) / rw));
			}
		}
	}
}

struct ClassicProblem{
	long long maximalValue(vector <int> cnt, vector <int> w, vector <int> v, int limit){
		int n = cnt.size();
		memset(dp, 0x3f, sizeof(dp));
		memset(dp1, 0, sizeof(dp1));
		dp[0] = 0;
		for(int i = 0; i < n; i++)
			arr[i] = pac(cnt[i], w[i], v[i]);
		sort(arr, arr + n);
		reverse(arr, arr + n);
		ll rsum = 0, ans = 0;
		for(int i = 0; i < n; i++) {
			rsum += arr[i].cnt * arr[i].w;
			ans += arr[i].cnt * arr[i].v;
			if(rsum > limit) {
				int rt = (rsum - limit - 1) / arr[i].w + 1;
				rsum -= rt * arr[i].w;
				ans -= rt * arr[i].v;
				vector <pac> arr1, arr2;
				for(int j = 0; j < i; j++)
					arr1.push_back(arr[j]);
				arr1.push_back(pac(arr[i].cnt - rt, arr[i].w, arr[i].v));
				arr2.push_back(pac(rt, arr[i].w, arr[i].v));
				for(int j = i + 1; j < n; j++)
					arr2.push_back(arr[j]);
				getdp(dp, 0, arr1);
				getdp(dp1, 1, arr2);
				for(int i = 1; i < maxv; i++)
					dp1[i] = max(dp1[i], dp1[i - 1]);
				ll rans = 0;
				for(int j = 0; j < maxv; j++)
					rans = max(rans, dp1[limit - rsum + j] - dp[j]);
				ans += rans;
				break;
			}

		}

		return ans ;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

