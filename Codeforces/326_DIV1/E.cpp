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

#ifndef N
#define N 200100
#endif
template <class t> struct segment_node{
	int be, en;
	t have, zhi[35], add, last, first;
};
template <class t> struct segment_tree{
	int l;
	segment_node <t> tree[N * 4];

	inline int gleft(int no) {return no << 1;}
	inline int gright(int no) {return (no << 1) + 1;}
	inline int gfa(int no) {return no >> 1;}
	inline segment_tree(){ l = 0; }

	void build(int no, int l, int r, t *a){
		if(l > r) r = l;
		if(l == r){
			tree[no].be = tree[no].en = l;
			tree[no].add = 0;
			tree[no].have = 0;
			tree[no].first = tree[no].last = a[l];
			if(a[l] == 0)
				tree[no].have = 0;
			return;
		}
		tree[no].be = l; tree[no].en = r;
		int mid = (l + r) / 2;
		build(gleft(no), l, mid, a);
		build(gright(no), mid + 1, r, a);
		pushup(tree + no, tree + gleft(no), tree + gright(no));
	}
	inline int gauss(int *arr, int l) {
		int ans = 0;
		for(int i = 0; i < 31 && ans < l; i++) {
			int no = ans;
			for(int j = ans; j < l; j++)
				if(arr[j] & (1 << i))
					no = j;
			if(! (arr[no] & (1 << i)))
				continue;
			swap(arr[no], arr[ans]);
			for(int j = 0; j < l; j++)
				if(j != ans && (arr[j] & (1 << i)))
					arr[j] ^= arr[ans];
			ans++;
		}
		return ans;
	}
	inline void pushup(segment_node <t> *no, segment_node <t> *le, segment_node <t> *ri) {
//		if(le->have == 31 || ri->have == 31) {
//			no->have = 31;
//			for(int i = 0; i < 31; i++)
//				no->zhi[i] = 1 << i;
//			return;
//		}
		int arr[70];
		for(int i = 0; i <le->have; i++)
			arr[i] = le->zhi[i];
		for(int j = 0; j < ri->have; j++)
			arr[j + le->have] = ri->zhi[j];
		arr[ri->have + le->have] = le->last ^ ri->first;
		no->have = gauss(arr, le->have + ri->have + 1);
		for(int i = 0; i < no->have; i++)
			no->zhi[i] = arr[i];
		no->last = ri->last;
		no->first = le->first;
	}
	inline void relax(segment_node <t> *no, segment_node <t> *le, segment_node <t> *ri) {
//		if(le->have < 31) {
//			for(int i = 0; i < max(le->have, 1); i++)
//				le->zhi[i] ^= no->add;
//			le->have = gauss(le->zhi, max(le->have, 1));
//		}
		le->last ^= no->add;
		le->add ^= no->add;
		le->first ^= no->add;
//		if(ri->have < 31) {
//			for(int i = 0; i < max(ri->have, 1); i++)
//				ri->zhi[i] ^= no->add;
//			gauss(ri->zhi, max(ri->have, 1));
//		}
		ri->last ^= no->add;
		ri->add ^= no->add;
		ri->first ^= no->add;
		no->add = 0;
	}
	void down(int l, int r, int no, t ranadd){
		if(l <= tree[no].be && r >= tree[no].en){
			tree[no].add ^= ranadd;
			tree[no].last ^= ranadd;
			tree[no].first ^= ranadd;
//			if(tree[no].have < 31) {
//				for(int i = 0; i < max(tree[no].have, 1); i++)
//					tree[no].zhi[i] ^= ranadd;
//				tree[no].have = gauss(tree[no].zhi, max(tree[no].have, 1));
//			}
			return;
		}
		if(tree[no].add && tree[no].be != tree[no].en)
			relax(tree + no, tree + gleft(no), tree + gright(no));
		int mid = (tree[no].be + tree[no].en) >> 1;
		if(r >= tree[no].be && l <= mid) down(l, r, gleft(no), ranadd);
		if(r >= mid + 1 && l <= tree[no].en) down(l, r, gright(no), ranadd);
		pushup(tree + no, tree + gleft(no), tree + gright(no));
	}
	segment_node <t> getsum(int l, int r, int no){
		if(l <= tree[no].be && r >= tree[no].en)
			return tree[no];
		if(tree[no].add && tree[no].be != tree[no].en)
			relax(tree + no, tree + gleft(no), tree + gright(no));
		segment_node <t> ans, le, ri;
		le.have = 0, ri.have = 0;
		le.first = ri.first = le.last = ri.last = 0;
		int mid = (tree[no].be + tree[no].en) >> 1;
		int flag = 0;
		if(r >= tree[no].be && l <= mid) {
			le = getsum(l, r, gleft(no));
			flag |= 1;
		}
		if(r >= mid + 1 && l <= tree[no].en) {
			ri = getsum(l, r, gright(no));
			flag |= 2;
		}
		if(flag == 3)
			pushup(&ans, &le, &ri);
		else if(flag == 2)
			ans = ri;
		else ans = le;
		return ans;
	}
};

const int maxn = N;
segment_tree <int> sgt;
int n, nq;
int arr[maxn];

int main() {

//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	scanf("%d%d", &n, &nq);
	for(int i = 0; i < n; i++)
		scanf("%d", arr + i);
	sgt.build(1, 0, n - 1, arr);
	for(int i = 0; i < nq; i++) {
		int cate;
		scanf("%d", &cate);
		if(cate == 1) {
			int l, r, num;
			scanf("%d%d%d", &l, &r, &num);
			l--, r--;
			sgt.down(l, r, 1, num);
		} else {
			int l, r;
			scanf("%d%d", &l, &r);
			l--; r--;
			segment_node <int> sgtnode = sgt.getsum(l, r, 1);
			sgtnode.zhi[sgtnode.have] = sgtnode.last;
			int ans = sgt.gauss(sgtnode.zhi, sgtnode.have + 1);
			printf("%I64d\n", (1ll << ans));
		}
	}

    return 0;
}
