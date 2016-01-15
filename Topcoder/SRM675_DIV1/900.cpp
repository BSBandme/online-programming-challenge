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
const double eps = 1e-10;
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

struct point {
	union {
		double co[2];
		struct {
			double x, y;
		};
	};
	double ang;
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
	double operator * (const point a) const {
		return x * a.x + y * a.y;
	}
	double operator % (const point a) {
		return x * a.y - y * a.x;
	}
	point operator * (double p) {
		point ans;
		ans.x = x * p;
		ans.y = y * p;
		return ans;
	}

};
double getang(double x, double y) {
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
	double ang;
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
	double sa = (lb.a - la.a) % (lb.a - la.b);
	double sb = (lb.b - la.b) % (lb.b - la.a);
	point ans = (lb.a * sb + lb.b * sa) *(1.0 / (sb + sa));
	return ans;
}

const int maxn = 55;
point orig[maxn];
double pr[maxn], pb[maxn];
int n;

double getport(point from, point to, point i) {
//	if(!(contain(line(from, to), i))) {
//		cerr << 1 << endl;
//	}
	if(to.x == from.x)
		return (i - from).y / (to - from).y;
	return (i - from).x / (to - from).x;
}
bool cmp(const pair <pdd, int> &a, const pair <pdd, int> &b) {
	return a.first.second < b.first.second;
}
double getl(double y1, double y2, double ratio) {
	return y2 * ratio + y1 * (1 - ratio);
}

struct BichromeSky2{
	double justdoit(int from, int to, double p[maxn]) {
		double ansl = 0.0;
		if(orig[from].x == orig[to].x)
			return 0.0;
		int must[maxn], cnt = 0;;
		memset(must, 0, sizeof(must));
		for(int i = 0; i < n; i++) if(i != from && i != to) {
			if(contain(line(orig[from], orig[to]), orig[i])) {
				double rp = getport(orig[from], orig[to], orig[i]);
				if(rp >= 0 && rp <= 1)
					continue;
				else
					must[i] = 1, cnt++;
			} else if(in(orig[from], orig[to], orig[i]))
				must[i] = 1, cnt++;
		}

//		printf("%d %d ", from, to);
//		for(int i = 0; i < n; i++)
//			printf("%d", must[i]);
//		puts("");

		if(cnt == 0)
			return 0;

		double prep = p[from] * p[to];
		pair <pdd,int> range[maxn];
		int lr = 0;
		for(int i = 0; i < n; i++) if(i != from && i != to){
			if(must[i]) {
				prep *= 1 - p[i];
				continue;
			}
			pdd k = mpr(1, 0);
			for(int j = 0; j < n; j++) if(must[j]){
				point interc = getinter(line(orig[j], orig[i]), line(orig[from], orig[to]));
				double rp = getport(orig[from], orig[to], interc);
				k.first = min(k.first, rp);
				k.second = max(k.second, rp);
			}
			k.first = max(0.0, k.first);
			k.second = min(1.0, k.second);
			range[lr++] = mpr(k, i);
		}

		sort(range, range + lr);
		double pl = 1.0;
		double origarea = (orig[to].x - orig[from].x) / 2;
		for(int i = 0; i < lr; i++) {
			pair <pdd,int> rk[maxn];
			int lrk = 0;
			for(int j = i + 1; j < lr; j++) if(jud(range[j].first.second, range[i].first.second) == 1) {
				rk[lrk++] = range[j];
			}
			sort(rk, rk + lrk, cmp);
			double pr = 1.0;
			for(int j = lrk - 1; j >= 0; j--) {
				double rp = pl * pr * (1 - p[range[i].second]) * (1 - p[rk[j].second]);
				ansl += prep * rp * (rk[j].first.second - range[i].first.first) * origarea * (getl(orig[from].y, orig[to].y, rk[j].first.second) + getl(orig[from].y, orig[to].y, range[i].first.first));
				pr *= p[rk[j].second];
			}
			ansl += prep * pl * pr * (1 - p[range[i].second]) * (range[i].first.second - range[i].first.first) * origarea * (getl(orig[from].y, orig[to].y, range[i].first.second) + getl(orig[from].y, orig[to].y, range[i].first.first));
			pl *= p[range[i].second];
		}
		return ansl;

	}

	double expectationOfArea(vector <int> x, vector <int> y, vector <int> prob){
//		reverse(x.begin(), x.end());
//		reverse(y.begin(), y.end());
//		reverse(prob.begin(), prob.end());
//		swap(x[1], x[4]);
//		swap(y[1], y[4]);
//		swap(prob[1], prob[4]);
		n = x.size();
		for(int i = 0; i < n; i++) {
			orig[i] = point(x[i], y[i]);
			pr[i] = prob[i] / 1000.0;
			pb[i] = 1 - pr[i];
		}

		double ans = 0.0;
		for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) if(i != j){
//			if(i == 0 && j == 5) {
//				cerr << 1 << endl;
//			}
			ans += justdoit(i, j, pr);
			ans += justdoit(i, j, pb);
		}

		return ans;
	}
	


};





// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

