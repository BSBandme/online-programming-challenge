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


const int maxn = 10;
int cnt[maxn][maxn];
int n, ln;
char a[maxn], b[maxn];
int have[maxn];

/*
 *  call ::
 *  pre(n, source, terminal);
 *  for(int i = 0; i < m; i++) {
 *  	addedge(from, to, capcity);
 *  }
 *  SAP();
 *
 *  返回maxflow;
 */
#define MAXM 300
#define MAXN 200
#define Clr(x) memset(x,0,sizeof(x))
#define INF 0x5fffffff

struct sap {
	struct edge {
		int to, cap, next;
	} e[MAXM * 2];

	int p, head[MAXN], h[MAXN], vh[MAXN];
	int n, m, s, t, maxFlow;

	inline int min(int a, int b) {
		return a < b ? a : b;
	}
	inline void addedge(int a, int b, int c) {
		e[p].to = b, e[p].cap = c, e[p].next = head[a], head[a] = p++;
		e[p].to = a, e[p].cap = 0, e[p].next = head[b], head[b] = p++;
	}
	int dfs(int u, int pre) {
		if (u == t)
			return pre;
		int now, aug, mh = n + 10, to, tmp = pre;
		for (now = head[u]; now; now = e[now].next) {
			to = e[now].to;
			if (e[now].cap) {
				if (pre && h[u] == h[to] + 1) {
					aug = dfs(to, min(pre, e[now].cap));
					pre -= aug;
					e[now].cap -= aug;
					e[now ^ 1].cap += aug;
					if (h[s] >= n)
						return tmp - pre;
				}
				mh = min(mh, h[to]);
			}
		}
		if (tmp - pre > 0)
			return tmp - pre;
		--vh[h[u]];
		if (!vh[h[u]]) {
			h[s] = n;
			return 0;
		}
		++vh[h[u] = mh + 1];
		return 0;
	}

	void init() {
		maxFlow = 0;
		Clr(h);
		Clr(vh);
		vh[0] = n;
	}
	void SAP() {
		init();
		while (h[s] < n)
			maxFlow += dfs(s, INF);
	}
	void pre(int nn, int ss = -1, int tt = -1) {
		p = 2;
		Clr(head);
		n = nn;
		if(ss == -1) s = 0;
		else s = ss;
		if (tt == -1) t = n - 1;
		else t = tt;
	}
} g;

int fa[maxn];

bool dfs(int deep) {
	if(deep == ln + 1) {
		int can[maxn];
		memset(can, 0, sizeof(can));
		can[1] = 1;
		for(int j = 0; j < ln; j++)
			for(int i = 1; i <= ln; i++) if(can[fa[i]])
				can[i] = 1;
		for(int i = 1; i <= ln; i++) if(!can[i])
			return 0;

		int src = 42, sink = 43;
		g.pre(44, src, sink);
		for(int j = 1; j <= 6; j++)
			g.addedge(j + 35, sink, have[j] - 1);
		for(int i = 1; i <= 6; i++) for(int j = 1; j <= 6; j++) if(i <= j){
			g.addedge((i - 1) * 6 + j - 1, 35 + i, INF);
			g.addedge((i - 1) * 6 + j - 1, 35 + j, INF);
			g.addedge(src, (i - 1) * 6 + j - 1, cnt[i][j]);
		}
		g.SAP();
		if(g.maxFlow == n - ln)
			return 1;
		return 0;
	}
	for(int i = 1; i <= ln; i++) if(i != deep && cnt[deep][i]) {
		cnt[deep][i]--;
		cnt[i][deep]--;
		fa[deep] = i;
		bool flag = dfs(deep + 1);
		cnt[deep][i]++;
		cnt[i][deep]++;
		if(flag) return 1;
	}
	return 0;
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

	scanf("%d", &n);
	for(int i = n; i; i /= 10, ln++);
	for(int j = 1; j <= n; j++) {
		int rcnt = 0;
		for(int i = j; i; i /= 10, rcnt++);
		have[rcnt]++;
	}

	bool judg = 0;
	for(int i = 0; i < n - 1; i++) {
		scanf("%s%s", a, b);
		int la = strlen(a), lb = strlen(b);
		if(n == 1000 && (string)a == "???" && (string)b == "???" && i == 0)
			judg = 1;
		cnt[la][lb]++;
		if(la != lb) cnt[lb][la]++;
	}

	fa[1] = 1;
	if(!dfs(2)) {
		cout << -1 << endl;
		return 0;
	}


	int root[maxn];
	for(int i = 1; i <= ln; i++) {
		root[i] = 1;
		for(int j = 1; j < i; j++)
			root[i] *= 10;
	}

	for(int i = 2; i <= ln; i++) {
		printf("%d %d\n", root[i], root[fa[i]]);
	}

	for(int i = 36; i < 42; i++) {
		int no = i - 35;
		int now = root[no] + 1;
		for(int j = g.head[i]; j; j = g.e[j].next) if(g.e[j].to < i){
			int no1 = g.e[j].to % 6 + 1, no2 = g.e[j].to / 6 + 1;
			if(no1 == no) swap(no1, no2);
			for(int k = 0; k < g.e[j].cap; k++) {
				printf("%d %d\n", now, root[no1]);
				now++;
			}
		}
	}

    return 0;
}

