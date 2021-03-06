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
#ifndef N
#define N 200100
#endif
template <class t> struct segment_node{
	int be, en;
	t multi, sum;
};
template <class t> struct segment_tree{
	int l;
	segment_node <t> tree[N * 4];

	inline int gleft(int no) {return no << 1;}
	inline int gright(int no) {return (no << 1) + 1;}
	inline int gfa(int no) {return no >> 1;}
	inline segment_tree(){ l = 0; }

	void build(int no, int l, int r, t orig = 0, t *a = NULL){
		if(l > r) r = l;
		if(l == r){
			tree[no].be = tree[no].en = l;
			tree[no].multi = 1;
			if(a == NULL) tree[no].sum = orig;
			else tree[no].sum = a[l];
			tree[no].multi = 1;
			return;
		}
		tree[no].be = l; tree[no].en = r;
		int mid = (l + r) / 2;
		build(gleft(no), l, mid, orig, a);
		build(gright(no), mid + 1, r, orig, a);
		tree[no].sum = tree[gleft(no)].sum + tree[gright(no)].sum;
		tree[no].sum %= mod;
		tree[no].multi = 1;
	}
	inline void relax(int no) {
		int le = gleft(no), ri = gright(no);
		tree[le].multi *= tree[no].multi; tree[le].multi %= mod;
		tree[le].sum *= tree[no].multi; tree[le].sum %= mod;
		tree[ri].multi *= tree[no].multi; tree[ri].multi %= mod;
		tree[ri].sum *= tree[no].multi; tree[ri].sum %= mod;
		tree[no].multi = 1;
	}
	void down(int l, int r, int no, t ranadd){
		if(l <= tree[no].be && r >= tree[no].en){
			tree[no].multi *= ranadd; tree[no].multi %= mod;
			tree[no].sum *= ranadd; tree[no].sum %= mod;
			return;
		}
		if(tree[no].multi != 1 && tree[no].be != tree[no].en) relax(no);
		int mid = (tree[no].be + tree[no].en) >> 1;
		if(r >= tree[no].be && l <= mid) down(l, r, gleft(no), ranadd);
		if(r >= mid + 1 && l <= tree[no].en) down(l, r, gright(no), ranadd);
		tree[no].sum = tree[gleft(no)].sum + tree[gright(no)].sum;
		tree[no].sum %= mod;
	}
	t getsum(int l, int r, int no){
		if(l > r) return 0;
		if(l <= tree[no].be && r >= tree[no].en)
			return tree[no].sum;
		if(tree[no].multi != 1 && tree[no].be != tree[no].en) relax(no);
		t ans = 0;
		int mid = (tree[no].be + tree[no].en) >> 1;
		if(r >= tree[no].be && l <= mid) ans += getsum(l, r, gleft(no));
		if(r >= mid + 1 && l <= tree[no].en) ans += getsum(l, r, gright(no));
		return ans % mod;
	}
	void set(int loc, int no, t ranadd){
		if(tree[no].be == tree[no].en) {
			tree[no].sum = ranadd;
			return;
		}
		if(tree[no].multi != 1 && tree[no].be != tree[no].en) relax(no);
		int mid = (tree[no].be + tree[no].en) >> 1;
		if(loc <= mid) set(loc, gleft(no), ranadd);
		else set(loc, gright(no), ranadd);
		tree[no].sum = tree[gleft(no)].sum + tree[gright(no)].sum;
		tree[no].sum %= mod;
	}
};

segment_tree <ll> sgt, xi;

const int maxn = 200100;
struct edge {
	int to, nxt;
} e[maxn * 2];
int head[maxn], le;
pii range[maxn];
int fa[maxn];
ll val[maxn], rv[maxn];
int prearr[maxn], lparr;
int nq, n;
int orig[maxn], du[maxn];
ll inv[maxn];

void dfs(int no) {
	prearr[lparr] = no;
	range[no].first = lparr++;
	for(int i = head[no]; i != -1; i = e[i].nxt) {
		dfs(e[i].to);
	}
	range[no].second = lparr - 1;
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

	memset(head, -1, sizeof(head));
	scanf("%I64d%d", val, &nq);
	n = 1;
	for(int i = 0; i < nq; i++) {
		int cate;
		scanf("%d", &cate);
		if(cate == 1) {
			int a, b;
			scanf("%d%d", &a, &b);

			e[le].to = n;
			e[le].nxt = head[a - 1];
			head[a - 1] = le++;

			val[n] = b;
			fa[n] = a - 1;
			orig[i] = n++;
		} else {
			int a;
			scanf("%d", &a);
			orig[i] = -(a - 1);
		}
	}
	dfs(0);
//	for(int i = 0; i < n; i++)
//		rv[range[i].first] = val[i];

	inv[1] = 1;
	for(int i = 2; i < maxn; i++)
		inv[i] = pow(i, mod - 2, mod);

	sgt.build(1, 0, n - 1, 0);
	xi.build(1, 0, n - 1, 1);

	sgt.set(range[0].first, 1, val[0]);
	for(int i = 0; i < nq; i++) {
		//
//		#ifdef DEBUG //......................................................................................................
//		for(int i = 0; i < n; i++)
//			printf("%I64d ", sgt.getsum(range[i].first, range[i].first, 1));
//		puts("");
//		for(int i = 0; i < n; i++)
//			printf("%I64d ", xi.getsum(range[i].first, range[i].first, 1));
//		puts("");
//		#endif //...........................................................................................................

		int no = orig[i];
		if(no > 0) {
			du[fa[no]]++;
			xi.down(range[fa[no]].first, range[fa[no]].second, 1, inv[du[fa[no]]] * (du[fa[no]] + 1) % mod);
			sgt.down(range[fa[no]].first, range[fa[no]].second, 1, inv[du[fa[no]]] * (du[fa[no]] + 1) % mod);
			sgt.set(range[no].first, 1, xi.getsum(range[no].first, range[no].first, 1) * val[no] % mod);
		} else {
			no = -no;
			ll ans = sgt.getsum(range[no].first, range[no].second, 1);
			ans *= pow(xi.getsum(range[no].first, range[no].first, 1), mod - 2, mod);
			ans %= mod;
			ans *= du[no] + 1;
			printf("%I64d\n",  ans % mod);
		}


	}

    return 0;
}
