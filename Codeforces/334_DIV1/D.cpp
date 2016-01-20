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
const long double eps = 1e-11;
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


struct point {
	union {
		long double co[2];
		struct {
			long double x, y;
		};
	};
	long double ang;
	point(double a = 0, double b = 0) : x(a), y(b) {ang = 0;}

	point operator + (const point a) {
		point ans;
		ans.x = x + a.x;
		ans.y = y + a.y;
		return ans;
	}
	point operator - (const point a) {
		point ans;
		ans.x = x - a.x;
		ans.y = y - a.y;
		return ans;
	}
	long double operator * (const point a) const {
		return x * a.x + y * a.y;
	}
	long double operator % (const point a) {
		return x * a.y - y * a.x;
	}
	point operator * (double p) {
		point ans;
		ans.x = x * p;
		ans.y = y * p;
		return ans;
	}

};
long double getang(long double x, long double y) {
	if(x == 0) {
		if(y > 0) return pi / 2;
		else return pi / 2 * 3;
	}
	if(x > 0) {
		if(y < 0) return atan(y / x) + pi * 2;
		else return atan(y / x);
	} else {
		return atan(y / x) + pi;
	}
}
double getang(point p) {
	return getang(p.x, p.y);
}
bool in(point cen, point a, point b) {
	return jud((a - cen) % (b - cen), 0) >= 0;
}

//判断两线段是否相交
bool inter(point a, point b, point c, point d){
	if ( min(a.x, b.x) > max(c.x, d.x) ||
	min(a.y, b.y) > max(c.y, d.y) ||
	min(c.x, d.x) > max(a.x, b.x) ||
	min(c.y, d.y) > max(a.y, b.y) ) return 0;
	double h, i, j, k;
	h = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	i = (b.x - a.x) * (d.y - a.y) - (b.y - a.y) * (d.x - a.x);
	j = (d.x - c.x) * (a.y - c.y) - (d.y - c.y) * (a.x - c.x);
	k = (d.x - c.x) * (b.y - c.y) - (d.y - c.y) * (b.x - c.x);
	return h * i <= eps && j * k <= eps;
}
struct line {
	point a, b;
	long double ang;
	line(point a1 = point(0, 0), point b1 = point(1, 0)): a(a1), b(b1) {
		point temp = b1 - a1;
		ang = getang(temp.x, temp.y);
	}
	friend bool contain(line l, point a);
	bool operator == (const line &rline) const {
		return contain(*this, rline.a) && contain(*this, rline.b);
	}
};
inline bool contain(line l, point a) {
	return jud((a - l.a) % (l.b - l.a), 0) == 0;
}
point getinter(line la, line lb) {
	long double sa = (lb.a - la.a) % (lb.a - la.b);
	long double sb = (lb.b - la.b) % (lb.b - la.a);
	point ans = (lb.a * sb + lb.b * sa) *(1.0 / (sb + sa));
	return ans;
}

const int maxn = 2010;
line orig[maxn];
int n;

line getline(int a, int b, int c) {
	line l;
	if(a == 0) {
		l.a.x = 0;
		l.a.y = 1.0 * c / b;
		l.b.x = 1;
		l.b.y = 1.0 * c / b;
	} else if(b == 0) {
		l.a.x = 1.0 * c / a;
		l.a.y = 0;
		l.b.x = 1.0 * c / a;
		l.b.y = 1;
	} else if(c == 0) {
		l.a.x = 1;
		l.a.y = -1.0 * a / b;
		l.b.x = -1;
		l.b.y = 1.0 * a / b;
	} else {
		l.a.x = 1.0 * c / a;
		l.a.y = 0;
		l.b.x = 0;
		l.b.y = 1.0 * c / b;
	}
	return l;
}

pair <long double, int> ri[maxn], le[maxn];
int cont[maxn];
int lle, lri;
ll ans;

int main() {

//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

	ios_base::sync_with_stdio(0);
	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	scanf("%d", &n);
	int rn = 0;
	for(int i = 0; i < n; i++) {
		int a, b, c;
		scanf("%d%d%d", &a, &b, &c);
		orig[i] = getline(a, b, c);
		if(c == 0) {
			rn++;
			cont[i] = 1;
		}
	}
	ans = rn * (rn - 1) / 2 * (n - rn);
	for(int i = 0; i < n; i++) if(!cont[i]){
		lle = 0;
		lri = 0;
		for(int j = 0; j < n; j++) if(j != i && !cont[j]) {
			point p = getinter(orig[j], orig[i]);
			point p1 = orig[j].a;
			if(jud(p.x, p1.x) == 0 && jud(p.y, p1.y) == 0)
				p1 = orig[j].b;
			long double f = (p1 - orig[i].b) % (orig[i].a - orig[i].b);
			f *= (point(0, 0) - orig[i].b) % (orig[i].a - orig[i].b);
			if(f * (p1 % p) > 0) {
				long double rang = -getang(p1 - p) + getang(point(0,0) - p);
				if(f > 0) rang = pi + rang;
				if(rang < 0) rang += 2 * pi;
				if(rang > 2 * pi) rang -= 2 * pi;
				if(jud(rang, 0) == 0 || jud(rang, pi) == 0)
					continue;
				ri[lri++] = mpr(rang, j);
			} else {
				long double rang = getang(p1 - p) - getang(point(0,0) - p);
				if(f > 0) rang = pi + rang;
				if(rang < 0) rang += 2 * pi;
				if(rang > 2 * pi) rang -= 2 * pi;
				if(jud(rang, 0) == 0 || jud(rang, pi) == 0)
					continue;

				le[lle++] = mpr(rang, j);
			}
		}

		sort(ri, ri + lri);
		sort(le, le + lle);
		for(int i = 0, j = lri - 1, rj = lri - 1; j >= 0 && i < lle; i++) {
			while(j >= 0 && jud(ri[j].first + le[i].first, pi) > 0)
				j--;
			while(rj >= 0 && jud(ri[rj].first + le[i].first, pi) >= 0)
				rj--;
			ans += j - rj;
		}
	}
	cout << ans << endl;

    return 0;
}
