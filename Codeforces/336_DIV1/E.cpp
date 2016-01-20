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
//template <typename t> inline int jud(t a, t b){
//	if(a < b) return -1;
//	if(a == b) return 0;
//	return 1;
//}

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

const int maxn = 100100;
pdd arr[maxn];
pdd cen;
int n, m;
pair <long double, int> inter[maxn];
int cnt[maxn], loc[maxn];
int bef[maxn], nxt[maxn];

pair <pair <double, double>, pair <double, double>> getinter(int no, double r) {
	double ra = (1 + arr[no].first * arr[no].first);
	double rb = 2 * (arr[no].second - cen.second) * arr[no].first - 2 * cen.first;
	double rc = (arr[no].second - cen.second) * (arr[no].second - cen.second) + cen.first * cen.first - r * r;
	if(rb * rb < 4 * ra * rc)
		return mpr(mpr(-1e100, -1e100), mpr(-1e100, -1e100));
	double delta = (rb * rb - 4 * ra * rc);
	delta = sqrt(delta);
	double x1 = (-rb + delta) / 2 / ra, x2 = (-rb - delta) / 2 / ra;
//	assert(jud(hypot(x1 - cen.first, arr[no].first * x1 + arr[no].second - cen.second), r) == 0);
//	assert(jud(hypot(x2 - cen.first, arr[no].first * x2 + arr[no].second - cen.second), r) == 0);
	return mpr(mpr(x1, arr[no].first * x1 + arr[no].second), mpr(x2, arr[no].first * x2 + arr[no].second));
}

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

#ifndef N
#define N 1 << 17
#endif
template <class t> struct tree_array{
	t num[N], n, bigwei;
	void upd(int no, t add){
		while(no <= n){
			num[no] += add;
			no += lowb(no);
		}
	}
	t que(int no){
		t ans = 0;
		while(no){
			ans += num[no];
			no -= lowb(no);
		}
		return ans;
	}
	 int getrank(t x){
		 int ans = 0, ranwei = bigwei; t rank = 0;
		 while(ranwei){
			 if(rank + num[ranwei + ans] < x){
				 ans += ranwei;
				 rank += num[ans];
			 }
			 ranwei >>= 1;
		 }
		 return ans + 1;
	 }
};

tree_array <int> ta;

int main() {


//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

	ios_base::sync_with_stdio(0);
	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	double rtime = 0;
	scanf("%d", &n);
	scanf("%lf%lf%d", &cen.first, &cen.second, &m);
	cen.first /= 1000, cen.second /= 1000;
	for(int i = 0; i < n; i++) {
		scanf("%lf%lf", &arr[i].first, &arr[i].second);
		arr[i].first /= 1000;
		arr[i].second /= 1000;
	}

	long double be = 0, en = 1e10;
	int lin;
	for(int ii = 0; ii < 100; ii++) {
		long double mid = (be + en) / 2;
		lin = 0;
		for(int i = 0; i < n; i++) {
			pair <pdd, pdd> temp = getinter(i, mid);
			if(temp.first.first < -1e50)
				continue;
			inter[lin++] = mpr(getang(temp.first.first - cen.first, temp.first.second - cen.second), i);
			inter[lin++] = mpr(getang(temp.second.first - cen.first, temp.second.second - cen.second), i);
		}
		sort(inter, inter + lin);
		memset(cnt, 0, sizeof(cnt));
		ta.n = lin + 10;
		for(int i = 0; i < lin + 10; i++)
			ta.num[i + 1] = 0;
		int rsum = 1;
		ll all = 0;
		for(int i = 0; i < lin; i++) {
			if(cnt[inter[i].second] == 0) {
				cnt[inter[i].second]++;
				loc[inter[i].second] = rsum;
				ta.upd(rsum, 1);
				rsum++;
			} else {
				cnt[inter[i].second]++;
				all += ta.que(rsum) - ta.que(loc[inter[i].second]);
				ta.upd(loc[inter[i].second], -1);
			}
		}
		if(all < m)
			be = mid;
		else
			en = mid;
	}

	lin = 0;
	for(int i = 0; i < n; i++) {
		pair <pdd, pdd> temp = getinter(i, be);
		if(temp.first.first < -1e50)
			continue;
		inter[lin++] = mpr(getang(temp.first.first - cen.first, temp.first.second - cen.second), i);
		inter[lin++] = mpr(getang(temp.second.first - cen.first, temp.second.second - cen.second), i);
	}
	sort(inter, inter + lin);
	memset(cnt, 0, sizeof(cnt));
	memset(bef, -1, sizeof(bef));
	memset(nxt, -1, sizeof(nxt));
	int head = -1;
	double ans = 0;
	int rcnt = 0;
	for(int i = 0; i < lin; i++) {
		int no = inter[i].second;
		if(cnt[no] == 0) {
			cnt[no]++;
			nxt[no] = head;
			if(head != -1)
				bef[head] = no;
			head = no;
		} else {
			double rx, ry;
			for(int j = head; j != no; j = nxt[j]) {
				rx = (arr[j].second - arr[no].second) / (arr[no].first - arr[j].first);
				ry = arr[no].first * rx + arr[no].second;
				ans += sqrt((cen.first - rx) * (cen.first - rx) + (cen.second - ry) * (cen.second - ry));
				rcnt++;
			}
			if(bef[no] != -1) nxt[bef[no]] = nxt[no];
			if(nxt[no] != -1) bef[nxt[no]] = bef[no];
			if(no == head) head = nxt[head];
		}
	}
	ans += (m - rcnt) * be;

	printf("%.10f\n", ans);


    return 0;
}


