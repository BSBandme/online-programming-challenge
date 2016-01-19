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

const int maxn = 10;
vll num;
int n;
ll mask = 0;
int cnt1[3000];
int gete[20][2];
int trans[1000];
mpii states;
int e[40000][40];
int arr[10];
map <pll, ll> mp[2];

void getStates(int deep, int already) {
	if(deep == n) {
		int temp = 0;
		for(int i = 0; i < n; i++)
			temp = temp * 10 + arr[i];
		int t = states.size();
		states[temp] = t;
		int rcnt = 0;
		for(int i = 0; i < n; i++) for(int j = i + 1; j < n; j++) {
			int tarr[10];
			memcpy(tarr, arr, sizeof(tarr));
			int p1 = tarr[i], p2 = tarr[j];
			if(p1 < p2) swap(p1, p2);
			e[t][rcnt] = 0;
			for(int k = 0; k < n; k++) {
				if(tarr[k] == p1)
					tarr[k] = p2;
				if(p1 != p2 && tarr[k] > p1)
					tarr[k]--;
				e[t][rcnt] = e[t][rcnt] * 10 + tarr[k];
			}
			e[t][rcnt] = states[e[t][rcnt]];
			rcnt++;
		}
	} else {
		if(arr[deep] != 0) {
			getStates(deep + 1, already);
			return;
		}
		int rmask = 0;
		for(int i = deep + 1; i < n; i++) if(arr[i] == 0)
			rmask |= 1 << i;
		already++;
		arr[deep] = already;
		for(int tmask = rmask; rmask >= 0; rmask = (rmask - 1) & tmask) {
			for(int t = 0; t < n; t++)
				if(rmask & (1 << t))
					arr[t] = already;
			getStates(deep + 1, already);
			for(int t = 0; t < n; t++)
				if(rmask & (1 << t))
					arr[t] = 0;
			if(rmask == 0) break;
		}
		already--;
		arr[deep] = 0;
		return;
	}
}

struct Gxor{

	int gauss() {
		int have = 0;
		for(int i = 0; i < (n - 1) * n / 2 && have < (int)num.size(); i++) {
			int no = have;
			for(int j = have; j < (int)num.size(); j++)
				if(num[j] & (1ll << i))
					no = j;
			if(!(num[no] & (1ll << i)))
				continue;
			mask |= 1ll << i;
			swap(num[have], num[no]);
			for(int j = 0; j < (int)num.size(); j++)
				if(j != have && (num[j] & (1ll << i)))
					num[j] ^= num[have];
			have++;
		}
		return have;
	}
	long long countsubs(vector <string> s){
		int p = s[0].size();
		mask = 0;
		for(n = 2; n < 10; n++)
			if(n * (n - 1) / 2 == p)
				break;
		for(int i = 0; i < 1 << n; i++)
			cnt1[i] = bitnum(i);
		memset(gete, 0, sizeof(gete));
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < i; j++)
				gete[i][0] += n - 1 - j;
			gete[i][1] = 1 << (n - 1 - i);
			gete[i][1]--;
		}
		for(int i = 0; i < 10; i++)
			trans[1 << i] = i;


		num.clear();
		for(int i = 0; i < (int)s.size(); i++) {
			ll temp = 0;
			for(int j = n * (n - 1) / 2 - 1; j >= 0; j--) {
				temp = temp * 2 + (s[i][j] == '1');
			}
			num.push_back(temp);
		}
		int have = gauss();
		ll ans = 0;
		if(have <= 22) {
			for(int i = 0; i < 1 << have; i++) {
				ll pmask = 0;
				for(int j = 0; j < have; j++) if(i & (1 << j))
					pmask ^= num[j];
				int all = 1, last = 0;
				for(; last != all;) {
					int newm = all ^ last;
					last = all;
					for(; newm; newm -= lowb(newm)) {
						int no = trans[lowb(newm)];
						all |= ((pmask >> gete[no][0]) & gete[no][1]) << (no + 1);
						for(int k = 0, yi = 0; k < no; k++, yi += n - k)
							all |= bool(pmask & (1ll << (yi + no - 1 - k))) << k;
					}

				}
				if(last == (1 << n) - 1)
					ans++;
			}
		} else {
			states.clear();
			getStates(0, 0);
			int now = 0, nxt = 1;
			ll tmask = 1ll << ((n - 1) * n / 2);
			tmask -= 1 + mask;
			mp[now].clear();
			mp[now][mpr(states.size() - 1, 0)] = 1;
			for(int i = 0, tcnt = 0; i < n * (n - 1) / 2; i++)
				if((1ll << i) & mask) {
					mp[nxt].clear();
					for(auto it: mp[now]) {
						ll no = it.first.first, rmask = it.first.second;
						ll rno = e[no][i], rrmask = rmask ^ (num[tcnt] & tmask);
						mp[nxt][mpr(rno, rrmask)] += it.second;
						mp[nxt][mpr(no, rmask)] += it.second;
					}
					swap(now, nxt);
					tcnt++;
				}
			for(auto it : mp[now]) {
				ll rmask = it.first.second, no = it.first.first;
				for(int i = 0; i < n * (n - 1) / 2; i++)
					if(rmask & (1ll << i))
						no = e[no][i];
				if(no == 0)
					ans += it.second;
			}
		}
		return ans << (num.size() - have);
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

