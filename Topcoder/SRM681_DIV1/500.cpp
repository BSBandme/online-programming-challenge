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

#define DEBUG
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

//const int maxn = 10000;
//const int yu = 4000;
//ll st[maxn][3], lst;
ll a, b;
//set <pair <int, pair <pii, int> > > s;
//set <pair <pii, int> > tst;

const int mod = MOD;
long long gnxt(ll num) {
	return	((num ^ a) + b) & ((1ll << 50) - 1);
}
long long gbef(ll num) {
	if(num < b) num += 1ll << 50;
	return (num - b) ^ a;
}
//void tadd(ll &a) {
//	a++;
//	if(a == maxn)
//		a = 0;
//}
//void tminus(ll &a) {
//	a--;
//	if(a < 0) a += maxn;
//}
//int bef(ll a) {
//	if(a == 0)
//		return maxn - 1;
//	return a - 1;
//}

struct LimitedMemorySeries2{
	int getSum(int n, long long x0, long long a, long long b){
		::a = a;
		::b = b;
		ll ans = 0;
		ll num = x0;
		for(ll i = 0; i < n; i++, num = gnxt(num)) {
			ll l = i - 1, r = i + 1, j = 1;
			ll numl = gbef(num), numr = gnxt(num);
			for(; l >= 0 && r < n && numl < num && numr < num; ) {
				j++; l--; r++;
				numl = gbef(numl);
				numr = gnxt(numr);
			}
			j--;
			ans += j ;
			ans %= mod;
		}
		return ans;
//		//
//		#ifdef DEBUG //......................................................................................................
//		ll rx = x0;
//		for(int i = 0; i < n; i++) {
//			cerr << rx << ' ';
//			rx = gnxt(rx);
//		}
//		cerr << endl;
//		#endif //...........................................................................................................
//		memset(st, 0, sizeof(st));
//
//		lst = 0;
//		st[lst][0] = -1;
//		st[lst][1] = 1ll << 60;
//		st[lst++][2] = -1;
//		tst.insert(mpr(mpr(0, x0), 0));
//		s.insert(mpr(1, mpr(mpr(0, x0), 0)));
//		tadd(rtst);
//		ll ans = 0;
//
//		for(int i = 1; i < n; i++) {
//			ll nxtnum = gnxt(tst.rbegin()->first.second);
//			while(tst.size() && tst.rbegin()->first.second < nxtnum) {
//				tst
//				ll rl = i - tst[bef(rtst)][0] - 1;
//				ans += min(rl, tst[bef(rtst)][2]);
//				tminus(rtst);
//			}
//			if(rtst != ltst && tst[bef(rtst)][1] == nxtnum) {
//				ans += min(tst[bef(rtst)][2], i - tst[bef(rtst)][0] - 1);
//				tst[bef(rtst)][2] = i - tst[bef(rtst)][0] - 1;
//				tst[bef(rtst)][0] = i;
//
//				while(ltst != rtst && i - tst[ltst][0] - 1 >= tst[ltst][2]) {
//					ans += tst[ltst][2];
//					if(tst[ltst][0] - st[lst - 1][0] >= yu) {
//						st[lst][0] = tst[ltst][0];
//						st[lst][1] = tst[ltst][1];
//						lst++;
//					}
//					tadd(ltst);
//				}
//				continue;
//			}
//			if(rtst == ltst) {
//				int range = tst[ltst][0];
//				for(; st[lst - 1][1] < nxtnum; ) {
//					range = st[lst - 1][0];
//					lst--;
//				}
//
////				if(st[lst - 1][2] == ) {
//					st[lst - 1][2] = 1;
//					for(int j = st[lst - 1][0] + 1; j < range; j++) {
//						ll tnxtnum = gnxt(st[lst - 1][1]);
//						if(j == 0) tnxtnum = x0;
//						while(lst && st[lst - 1][1] <= tnxtnum)
//							lst--;
//						st[lst][0] = j;
//						st[lst][1] = tnxtnum;
//						st[lst++][2] = 1;
//					}
//					while(st[lst - 1][1] < nxtnum)
//						lst--;
////					st[lst - 1][2] = 0;
////				}
//				tst[rtst][2] = i - st[lst - 1][0] - 1;
//
//			} else tst[rtst][2] = i - tst[bef(rtst)][0] - 1;
//			tst[rtst][0] = i;
//			tst[rtst][1] = nxtnum;
//			tadd(rtst);
//			while(ltst != rtst && i - tst[ltst][0] - 1 >= tst[ltst][2]) {
//				if(ltst == 9998) {
//					ltst = 9998;
//				}
//				ans += tst[ltst][2];
//				if(tst[ltst][0] - st[lst - 1][0] >= yu) {
//					st[lst][0] = tst[ltst][0];
//					st[lst][1] = tst[ltst][1];
//					lst++;
//				}
//				tadd(ltst);
//			}
//		}
//
//		tst[rtst][0] = n;
//		for(int i = ltst; i != rtst; i++) {
//			ans += min(tst[i][2], n - tst[i][0] - 1);
//		}
//		return ans % MOD ;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

