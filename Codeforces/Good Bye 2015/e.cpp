// #define DEBUG
//#define USEPB_DS
#define USETR1
#define CPPELEVEN
#define GPP

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

using namespace std;

#ifndef CPPELEVEN
#ifdef USETR1
#include <tr1/unordered_map>
#include <tr1/unordered_set>
using namespace tr1;
#endif
#else
#include <unordered_map>
#include <unordered_set>
#endif

#ifdef USEPB_DS
#include <ext/pb_ds/priority_queue.hpp>
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
// binomial_heap_tag, rc_binomial_heap_tag, thin_heap_tag, binary_heap_tag
typedef __gnu_pbds::priority_queue<int, greater<int>, pairing_heap_tag> pq_type;
// splay_tree_tag, ov_tree_tag
typedef tree <int, null_type, less <int>, rb_tree_tag, tree_order_statistics_node_update> tree_type;
#endif

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
const long double eps = 1e-10;
const int step[8][2] = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}, {-1, 1}, {1, 1}, {1, -1}, {-1, -1}};

template <class T> inline T abs1(T a) {return a < 0 ? -a : a;}

#ifndef CPPELEVEN
template <class T> inline T max1(T a, T b) { return b < a ? a : b; }
template <class T> inline T max1(T a, T b, T c) { return max1(max1(a, b), c); }
template <class T> inline T max1(T a, T b, T c, T d) { return max1(max1(a, b, c), d); }
template <class T> inline T max1(T a, T b, T c, T d, T e) { return max1(max1(a, b, c, d), e); }
template <class T> inline T min1(T a, T b) { return a < b ? a : b; }
template <class T> inline T min1(T a, T b, T c) { return min1(min1(a, b), c); }
template <class T> inline T min1(T a, T b, T c, T d) { return min1(min1(a, b, c), d); }
template <class T> inline T min1(T a, T b, T c, T d, T e) { return min1(min1(a, b, c, d), e); }
#else
template <typename t, typename t1>
t min1(t a, t1 b) { return a < b ? a : b; }
template <typename t, typename... arg>
t min1(t a, arg... arr) { return min1(a, min1(arr...)); }
template <typename t, typename t1>
t max1(t a, t1 b) { return a > b ? a : b; }
template <typename t, typename... arg>
t max1(t a, arg... arr) { return max1(a, max1(arr...)); }
#endif

inline int jud(double a, double b){
	if(abs(a) < eps && abs(b) < eps) return 0;
	else if(abs1(a - b) / max(abs1(a), abs1(b)) < eps) return 0;
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
	if(na == 0) return 0;
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
#ifdef GPP
inline int bitnum(ui nValue) { return __builtin_popcount(nValue); }
inline int bitnum(int nValue) { return __builtin_popcount(nValue); }
inline int bitnum(ull nValue) { return __builtin_popcount(nValue) + __builtin_popcount(nValue >> 32); }
inline int bitnum(ll nValue) { return __builtin_popcount(nValue) + __builtin_popcount(nValue >> 32); }
inline int bitmaxl(ui a) { if(a == 0) return 0; return 32 - __builtin_clz(a); }
inline int bitmaxl(int a) { if(a == 0) return 0; return 32 - __builtin_clz(a); }
inline int bitmaxl(ull a) { int temp = a >> 32; if(temp) return 32 - __builtin_clz(temp) + 32; return bitmaxl(int(a)); }
inline int bitmaxl(ll a) { int temp = a >> 32; if(temp) return 32 - __builtin_clz(temp) + 32; return bitmaxl(int(a)); }
#else
#endif

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

//....................密..........封..........线..........下..........禁..........止..........hack...............................................

const int maxn = 200100;
int arr[maxn];
int a, b, c;
multiset <int> ms;
int n;
int ans;

int main() {

//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

//	ios_base::sync_with_stdio(0);
	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	scanf("%d", &n);
	cin >> a >> b >> c;
	if(a > b) swap(a, b);
	if(a > c) swap(a, c);
	if(b > c) swap(b, c);
	for(int i = 0; i < n; i++) {
		scanf("%d", arr + i);
		ms.insert(arr[i]);
	}
	if(*ms.rbegin() > a + b + c) {
		cout << -1 << endl;
		return 0;
	}
	while(ms.size() && *ms.rbegin() > b + c) {
		auto it = ms.end();
		it--;
		ms.erase(it);
		ans++;
	}
	while(ms.size() && *ms.rbegin() > a + c) {
		auto it = ms.end();
		it--;
		ms.erase(it);
		ans++;
		it = ms.upper_bound(a);
		if(it != ms.begin()) {
			it--;
			ms.erase(it);
		}
	}
	while(ms.size() && *ms.rbegin() > max(c, a + b)) {
		auto it = ms.end();
		it--;
		ms.erase(it);
		ans++;
		it = ms.upper_bound(b);
		if(it != ms.begin()) {
			it--;
			ms.erase(it);
		}
	}
	if(a + b > c) {
		while(ms.size() && *ms.rbegin() > c) {
			auto it = ms.end();
			it--;
			ms.erase(it);
			ans++;
			it = ms.upper_bound(c);
			if(it != ms.begin()) {
				it--;
				ms.erase(it);
			}
		}
	}
	int flag = 0;
	for (; ms.size(); ) {
		flag = 0;
		auto it = ms.upper_bound(c);
		it--;
		ms.erase(it);
		if(ms.size()) {
			it = ms.upper_bound(b);
			if(it == ms.begin()) {
				it = ms.upper_bound(a + b);
			}
			if(it != ms.begin()) {
				it--;
				ms.erase(it);
			}
			if(ms.size()) {
				it = ms.upper_bound(a);
				if(it != ms.begin()) {
					it--;
					ms.erase(it);
				}
			}
		}
		ans++;
	}
//	int flag = 0;
//	for(; ms.size(); ) {
//		flag = 0;
//		auto ita = ms.upper_bound(a);
//		auto itb = ms.upper_bound(b);
//		auto itc = ms.upper_bound(c);
//		if(ita != ms.begin()) flag |= 1;
//		if(itb != ita) flag |= 2;
//		if(itc != itb) flag |= 4;
//		if(flag == 7) {
//			if(flag & 1) {
//				ita--;
//				ms.erase(ita);
//			}
//			if(flag & 2) {
//				itb--;
//				ms.erase(itb);
//			}
//			if(flag & 4) {
//				itc--;
//				ms.erase(itc);
//			}
//			ans++;
//		} else if(bitnum(flag) == 1){
//			int ra, rb;
//			if(flag == 1) {
//				ra = a;
//				rb = b + c;
//			} else if(flag == 2) {
//				ra = b;
//				rb = a + c;
//			} else {
//				ra = a + b;
//				rb = c;
//			}
//			if(ra > rb) swap(ra, rb);
//			ita = ms.upper_bound(ra);
//			itb = ms.upper_bound(rb);
//			if(itb == ita || ita == ms.begin()) {
//				int cnt = 0, allcnt = 0;
//				ita = ms.upper_bound(a);
//				if(ita != ms.begin()) cnt++;
//				ita = ms.upper_bound(b);
//				if(ita != ms.begin()) cnt++;
//				ita = ms.upper_bound(c);
//				if(ita != ms.begin()) cnt++;
//				if(cnt < 2) break;
//				ita--;
//				while(ita != ms.begin()) {
//					itb = ita;
//					ita--;
//					ms.erase(itb);
//					allcnt++;
//				}
//				ms.erase(ita);
//				allcnt++;
//				if(allcnt) ans += (allcnt - 1) / cnt + 1;
//				break;
//			} else {
//				ita--; itb--;
//				ms.erase(ita);
//				ms.erase(itb);
//				ans++;
//			}
//		} else if(bitnum(flag) == 2){
//			int ra1, rb1, ra2, rb2;
//			if(flag == 3) {
//				ra1 = a, rb1 = b + c;
//				ra2 = a + c, rb2 = b;
//			} else if(flag == 5) {
//				ra1 = a, rb1 = b + c;
//				ra2 = a + b, rb2 = c;
//			} else {
//				ra1 = b, rb1 = a + c;
//				ra2 = a + b, rb2 = c;
//			}
//			if(ra1 > rb1) swap(ra1, rb1);
//			if(ra2 > rb2) swap(ra2, rb2);
//			auto ita = ms.upper_bound(ra1);
//			auto itb = ms.upper_bound(rb1);
//			if(ita == itb || ita == ms.begin()) {
//				ita = ms.upper_bound(ra2);
//				itb = ms.upper_bound(rb2);
//			}
//			assert(ita != itb && ita != ms.begin());
//			ita--;
//			itb--;
//			ms.erase(ita);
//			ms.erase(itb);
//			ans++;
//		}
//	}

//	if(bitnum(flag) == 2) {
//	for(; ms.size(); ) {
//		auto ita = ms.upper_bound(a);
//		auto itb = ms.upper_bound(b);
//		auto itc = ms.upper_bound(c);
//		if(ita == itb && itb == itc && ita == ms.begin()) {
//			flag = 3;
//			break;
//		}
//		if(ita == itb && itb == itc) {
//			flag = 0;
//			break;
//		}
//		if(ita == ms.begin() && itb == itc) {
//			flag = 1;
//			break;
//		}
//		if(ita == ms.begin() && itb == ms.begin()) {
//			flag = 2;
//			break;
//		}
//		if(ita == ms.begin()) {
//			ms.erase(itc);
//			ms.erase(itb);
//		} else {
//			ms.erase(ita);
//			if(itc == itb)
//				ms.erase(itb);
//			else
//				ms.erase(itc);
//		}
//		ans++;
//	}
//	}
//
//
//	if(flag == 3) {
//		ans += ms.size();
//	} else {
//		if(flag == 0) b += c;
//		else if(flag == 1) {
//			a += c;
//			swap(a, b);
//		} else {
//			a += b;
//			b = c;
//			if(a > b) swap(a, b);
//		}
//		for(; ms.size(); ) {
//			auto ita = ms.upper_bound(a);
//			auto itb = ms.upper_bound(b);
//			if(ita == ms.begin() && itb == ita) break;
//			ita--; itb--;
//			ms.erase(ita);
//			ms.erase(itb);
//			ans++;
//		}
//		ans += ms.size();
//	}

	cout << ans + ms.size() << endl;

    return 0;
}
