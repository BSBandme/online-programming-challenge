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

const int mod = MOD;
const int maxn = 1010;
int arr[maxn];
int maxv[maxn];
int c[110][110];
int temp[400], temp1[400];

vi getdp(int le, int ri, vi already) {
	if(le >= ri - 1)
		return vi(1, 1);
	vi alle, alri;
	int mid = (le + ri) >> 1;
	int en = 110;
	int cntl = 0, cntr = 0;
	for(int i = 0; i < (int)already.size(); i++)  {
		if(already[i] < mid) {
			alle.push_back(already[i]);
			cntl++;
			if(cntl == mid - le) en = min(en, i);
		}
		else {
			alri.push_back(already[i]);
			cntr++;
			if(cntr == ri - mid) en = min(en, i);
		}
	}
	int alll = mid - le - cntl, allr = ri - mid - cntr;
	vi ansle = getdp(le, mid, alle);
	vi ansri = getdp(mid, ri, alri);
	vi ans(maxv[ri - le] + 10, 0);
	memset(temp1, 0, sizeof(temp1));
	if(alll == 0 || allr == 0) {
		temp1[en + 1] = 1;
	} else for(int i = min(alll, allr); i < alll + allr; i++) {
		add(temp1[i + cntl + cntr], c[i - 1][alll - 1]);
		add(temp1[i + cntl + cntr], c[i - 1][allr - 1]);
	}
	for(int i = 0; i < (int)ansle.size() + ansri.size(); i++)
		temp[i] = 0;
	for(int i = 0; i < (int)ansle.size(); i++)
		for(int j = 0; j < (int)ansri.size(); j++)
			add(temp[i + j], 1ll * ansle[i] * ansri[j] % mod);

	if(alll && allr) {
		for(int i = 0; i < (int)ansle.size() + ansri.size(); i++) if(temp[i]) {
			for(int j = min(alll, allr); j < alll + allr; j++) {
				add(ans[j + i + cntl + cntr], 1ll * temp1[j + cntl + cntr] * temp[i] % mod);
			}
		}
	} else {
		for(int i = 0; i < (int)ansle.size() + ansri.size(); i++) if(temp[i]){
			ans[en + 1 + i] = temp[i];
		}
	}

	return ans;
}

struct BearSorts{
	int index(vector <int> arr){
		for(int i = 0; i < 110; i++) {
			c[i][0] = 1;
			for(int j = 1; j <= i; j++) {
				c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
				c[i][j] %= mod;
			}
		}
		int n = arr.size();
		for(int i =0 ; i < n; i++)
			arr[i]--;
		maxv[1] = 1;
		for(int i = 2; i < maxn; i++) {
			maxv[i] = maxv[i / 2] + maxv[i - i / 2] + i - 1;
		}
		vi ans;
		ans = getdp(0, n, arr);
		int is = 0;
		for(int j = 0; j < (int)ans.size(); j++)
			if(ans[j]) is = j;
		ans = getdp(0, n, vi());
		ll sum = 0;
		for(int i = 0; i < is; i++)
			add(sum, ans[i]);
		vi tarr;
		int used[maxn];
		memset(used, 0, sizeof(used));
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < arr[i]; j++) if(used[j] == 0) {
				tarr.push_back(j);
				ans = getdp(0, n, tarr);
				add(sum, ans[is]);
				tarr.pop_back();
			}
			tarr.push_back(arr[i]);
			used[arr[i]] = 1;
		}
		return (sum + 1) % mod;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

