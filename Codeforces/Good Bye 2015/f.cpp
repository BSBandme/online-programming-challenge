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

const int mod = MOD;
const int maxn = 	500100;
char str[maxn];
int n, m, l, ln, rn, lm, rm, pln, prn, plm, prm;
int dx, dy, maxdx, maxdy, mindx, mindy;
int rmaxdx, rmaxdy, rmindx, rmindy, rdx, rdy;
ll ans;

void move(int no) {
	if(str[no] == 'L') dx--;
	else if(str[no] == 'R') dx++;
	else if(str[no] == 'U') dy--;
	else dy++;
	if(dx > maxdx) {
		rm--;
		ans += 1ll * (no + 1) * (rn - ln + 1) % mod;
	}
	if(dx < mindx) {
		lm++;
		ans += 1ll * (no + 1) * (rn - ln + 1) % mod;
	}
	if(dy > maxdy) {
		rn--;
		ans += 1ll * (no + 1) * (rm - lm + 1) % mod;
	}
	if(dy < mindy) {
		ln++;
		ans += 1ll * (no + 1) * (rm - lm + 1) % mod;
	}
	ans %= mod;
	maxdx = max(maxdx, dx);
	mindx = min(mindx, dx);
	maxdy = max(maxdy, dy);
	mindy = min(mindy, dy);
}

inline ll sum2(ll be, ll en) {
	be--;
	ll ans = en * (en + 1) * (2 * en + 1) / 6 % mod;
	if(be > 0)
		ans -= be * (be + 1) * (2 * be + 1) / 6 % mod;
	ans += mod;
	return ans % mod;
}
inline ll sum1(ll be, ll en) {
	be--;
	ll ans = (en + 1) * en / 2 % mod;
	if(be)
		ans -= be * (be + 1) / 2 % mod;
	ans += mod;
	return ans % mod;
}

ll getkuai(ll from, ll to, ll x, ll ln, ll rn, ll lm, ll rm, ll dx, ll dy, ll no){
	if(dx == 0)
		return 0;
	if((rm - x) / dx < 0 || (lm - x) / dx < 0)
		return 0;
	ll be = 0, en = 0;
	if(dx > 0) {
		be = (lm - x - 1) / dx + 1;
		en = (rm - x) / dx;
	} else {
		be = (rm + 1 - x) / dx + 1;
		en = (lm - x) / dx;
	}
	if(be != 1 || be > en)
		return 0;
	if(to + dy > rn) to -= (to + dy - rn);
	if(from + dy < ln) from += (ln - from - dy);
	if(from > to) return 0;
	ll rans = 1ll *  (en - be + 1) * (((be + en) * l % mod + 2 * (no + 1)) % mod) % mod * (to - from + 1) % mod;
	if(rans % 2) rans = (rans + MOD) / 2;
	else rans /= 2;
	ll beto = be * dy + to - rn, ento = en * dy + to - rn;
	if(beto > 0 || ento > 0) {
		if(beto < 0) beto -= (beto / abs(dy) - 1) * abs(dy);
		if(ento < 0) ento -= (ento / abs(dy) - 1) * abs(dy);
		beto += rn;
		ento += rn;
		ll rbe, ren;
		if(dy == 0) {
			rbe = be, ren = en;
		} else {
			rbe = (beto - to) / dy, ren = (ento - to) / dy;
		}
		rans -= (
					sum2(rbe, ren) * l % mod * dy % mod +
					(dy * (no + 1) + (to - rn) * l) % mod * sum1(rbe, ren) % mod +
					(ren - rbe + 1) * (no + 1) % mod * (to - rn) % mod
				) % mod;
	}
	beto = be * dy + (from - 1) - rn, ento = en * dy + (from - 1) - rn;
	if(beto > 0 || ento > 0) {
		if(beto < 0) beto -= (beto / abs(dy) - 1) * abs(dy);
		if(ento < 0) ento -= (ento / abs(dy) - 1) * abs(dy);
		beto += rn;
		ento += rn;
		ll rbe, ren;
		if(dy == 0) {
			rbe = be, ren = en;
		} else {
			rbe = (beto - (from - 1)) / dy, ren = (ento - (from - 1)) / dy;
		}
		rans += (
					sum2(rbe, ren) * l % mod * dy % mod +
					(dy * (no + 1) + ((from - 1) - rn) * l) % mod * sum1(rbe, ren) % mod +
					(ren - rbe + 1) * (no + 1) % mod * ((from - 1) - rn) % mod
				) % mod;
	}
	beto = ln - (be * dy + from), ento = ln - (en * dy + from);
	if(beto > 0 || ento > 0) {
		if(beto < 0) beto -= (beto / abs(dy) - 1) * abs(dy);
		if(ento < 0) ento -= (ento / abs(dy) - 1) * abs(dy);
		beto = ln - beto;
		ento = ln - ento;
		ll rbe, ren;
		if(dy == 0) {
			rbe = be, ren = en;
		} else {
			rbe = (beto - from) / dy, ren = (ento - from) / dy;
		}
		rans -= (
					-sum2(rbe, ren) * l % mod * dy % mod +
					((ln - from) * l - dy * (no + 1)) % mod * sum1(rbe, ren) % mod +
					(ren - rbe + 1) * (no + 1) % mod * (ln - from) % mod
				) % mod;
	}
	beto = ln - (be * dy + (to + 1)), ento = ln - (en * dy + (to + 1));
	if(beto > 0 || ento > 0) {
		if(beto < 0) beto -= (beto / abs(dy) - 1) * abs(dy);
		if(ento < 0) ento -= (ento / abs(dy) - 1) * abs(dy);
		beto = ln - beto;
		ento = ln - ento;
		ll rbe, ren;
		if(dy == 0) {
			rbe = be, ren = en;
		} else {
			rbe = (beto - (to + 1)) / dy, ren = (ento - (to + 1)) / dy;
		}
		rans += (
					-sum2(rbe, ren) * l % mod * dy % mod +
					((ln - (to + 1)) * l - dy * (no + 1)) % mod * sum1(rbe, ren) % mod +
					(ren - rbe + 1) * (no + 1) % mod * (ln - (to + 1)) % mod
				) % mod;
	}
	rans %= mod;
	if(rans < 0) rans += mod;

	return rans;
}

void move1(int no) {
	if(str[no] == 'L') rdx--;
	else if(str[no] == 'R') rdx++;
	else if(str[no] == 'U') rdy--;
	else rdy++;
	if(rdx > rmaxdx) {
		ans += getkuai(pln, prn, prm, ln, rn, lm, rm, -dx, -dy, no) % mod;
		prm--;
	}
	if(rdx < rmindx) {
		ans += getkuai(pln, prn, plm, ln, rn, lm, rm, -dx, -dy, no) % mod;
		plm++;
	}
	if(rdy > rmaxdy) {
		ans += getkuai(plm, prm, prn, lm, rm, ln, rn, -dy, -dx, no) % mod;
		prn--;
	}
	if(rdy < rmindy) {
		ans += getkuai(plm, prm, pln, lm, rm, ln, rn, -dy, -dx, no) % mod;
		pln++;
	}
	ans %= mod;
	rmaxdx = max(rmaxdx, rdx);
	rmindx = min(rmindx, rdx);
	rmaxdy = max(rmaxdy, rdy);
	rmindy = min(rmindy, rdy);
}


int main() {

//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

	ios_base::sync_with_stdio(0);
	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	scanf("%d%d%d", &l, &n, &m);
	pln = ln = 0, prn = rn = n - 1;
	plm = lm = 0, prm = rm = m - 1;
	scanf("%s", str);
	for(int i = 0; i < l; i++) {
		move(i);
		if(rm < lm || rn < ln) break;
	}
	if(rm < lm || rn < ln) {
		cout << ans << endl;
		return 0;
	}
	if(dx == 0 && dy == 0) {
		cout << -1 << endl;
		return 0;
	}
	for(int i = 0; i < l; i++) {
		move1(i);
	}
	cout << ans % mod << endl;

    return 0;
}
