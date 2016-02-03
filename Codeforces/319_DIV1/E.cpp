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

#ifdef DEBUG
const int maxn = 50100;
#else
const int maxn = 500010;
#endif
pii e[maxn];
vi eq[maxn];
int n, m, k, nq;
pii que[maxn];
pair <int*, int> store[maxn * 50];
int fa[50][maxn], sz[50][maxn];
int deep[50][maxn];
int ecol[maxn];
int ls;

inline void save(int &a) {
	store[ls++] = mpr(&a, a);
}
pii getfa(int no, int col) {
	if(no == fa[col][no])
		return mpr(no, 0);
	pii rans = getfa(fa[col][no], col);
	save(fa[col][no]);
	fa[col][no] = rans.first;
	save(deep[col][no]);
	deep[col][no] ^= rans.second;
	return mpr(fa[col][no], deep[col][no]);
}
int uni(int a, int b, int col) {
	pii faa = getfa(a, col), fab = getfa(b, col);
	if(faa.first == fab.first) {
		if(faa.second == fab.second)
			return -1;
		return 0;
	}
	if(sz[col][faa.first] < sz[col][fab.first]) {
		swap(faa, fab);
		swap(a, b);
	}
	save(sz[col][faa.first]);
	sz[col][faa.first] += sz[col][fab.first];
	save(fa[col][fab.first]);
	fa[col][fab.first] = faa.first;
	if(fab.second == faa.second) {
		save(deep[col][fab.first]);
		deep[col][fab.first] = 1;
	}
	return 1;
}


#ifndef N
#define N maxn
#endif
struct segment_node{
	int be, en;
	vector <pii> need;
};
struct segment_tree{
	int l;
	segment_node tree[N * 4];

	inline int gleft(int no) {return no << 1;}
	inline int gright(int no) {return (no << 1) + 1;}
	inline int gfa(int no) {return no >> 1;}
	inline segment_tree(){ l = 0; }

	void build(int no, int l, int r){
		if(l > r) r = l;
		if(l == r){
			tree[no].be = tree[no].en = l;
			return;
		}
		tree[no].be = l; tree[no].en = r;
		int mid = (l + r) / 2;
		build(gleft(no), l, mid);
		build(gright(no), mid + 1, r);
	}
	void down(int l, int r, int no, pii ranadd){
		if(l > r) return;
		if(l <= tree[no].be && r >= tree[no].en){
			tree[no].need.push_back(ranadd);
			return;
		}
		int mid = (tree[no].be + tree[no].en) >> 1;
		if(r >= tree[no].be && l <= mid) down(l, r, gleft(no), ranadd);
		if(r >= mid + 1 && l <= tree[no].en) down(l, r, gright(no), ranadd);
	}
	void getans(int no) {
		int rls = ls;
		for(int i = 0; i < (int)tree[no].need.size(); i++) {
			int eno = tree[no].need[i].first;
			int a = e[eno].first, b = e[eno].second, col = tree[no].need[i].second;
			uni(a, b, col);
		}
		tree[no].need.clear();
		if(tree[no].be == tree[no].en) {
			int lno = tree[no].be, eno = que[lno].first, col = que[lno].second;
			int loc = find(lno, eq[eno].begin(), eq[eno].size());
			int nxt = nq;
			if(loc != (int)eq[eno].size() - 1)
				nxt = eq[eno][loc + 1];
			if(uni(e[eno].first, e[eno].second, col) != -1) {
				ecol[eno] = col;
				puts("YES");
			} else
				puts("NO");
			if(ecol[eno] != -1)
				down(lno + 1, nxt - 1, 1, mpr(eno, ecol[eno]));

			//
			#ifdef DEBUG //......................................................................................................
			cerr << "//col" << endl;
			for(int i = 0; i < k; i++) {
				for(int j = 0; j < n; j++) {
					cerr << getfa(j, i).second << ' ';
				}
				cerr << endl;
			}
			cerr << endl;
			#endif //...........................................................................................................

		} else {
			getans(gleft(no));
			getans(gright(no));
		}

		while(ls > rls) {
			ls--;
			*store[ls].first = store[ls].second;
		}
	}
};

segment_tree sgt;

int main() {

//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

	ios_base::sync_with_stdio(0);
	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
//	freopen("output.txt", "w", stdout);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	scanf("%d%d%d%d", &n, &m, &k, &nq);
	for(int i = 0; i < m; i++) {
		int a, b;
		scanf("%d%d", &a, &b);
		a--; b--;
		e[i] = mpr(a, b);
	}
	for(int i = 0; i < nq; i++) {
		scanf("%d%d", &que[i].first, &que[i].second	);
		que[i].first--;
		que[i].second--;
		eq[que[i].first].push_back(i);
		ecol[i] = -1;
	}
	for(int i = 0; i < k; i++)
		for(int j = 0; j < n; j++) {
			fa[i][j] = j;
			sz[i][j] = 1;
		}
	sgt.build(1, 0, nq - 1);
	sgt.getans(1);
	// adadf23

    return 0;
}
