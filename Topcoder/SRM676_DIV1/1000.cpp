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
#define MAXM 10000
#define MAXN 200
#define Clr(x) memset(x,0,sizeof(x))
#define INF 0x3ffffffffffffll

struct sap {
	struct edge {
		ll to, cap, next;
	} e[MAXM * 2];

	ll p, head[MAXN], h[MAXN], vh[MAXN];
	ll n, m, s, t, maxFlow;

	inline ll min(ll a, ll b) {
		return a < b ? a : b;
	}
	inline void addedge(ll a, ll b, ll c) {
		e[p].to = b, e[p].cap = c, e[p].next = head[a], head[a] = p++;
		e[p].to = a, e[p].cap = 0, e[p].next = head[b], head[b] = p++;
	}
	ll dfs(ll u, ll pre) {
		if (u == t)
			return pre;
		ll now, aug, mh = n + 10, to, tmp = pre;
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
	void pre(ll nn, ll ss = -1, ll tt = -1) {
		p = 2;
		Clr(head);
		n = nn;
		if(ss == -1) s = 0;
		else s = ss;
		if (tt == -1) t = n - 1;
		else t = tt;
	}
	vi getcut() {
		vi ans;
		int q[220], lq = 0;
		int cont[220];
		memset(cont, 0, sizeof(cont));
		q[lq++] = s;
		cont[s] = 1;
		for(int i = 0; i < lq; i++) {
			int no = q[i];
			for(int j = head[no]; j; j = e[j].next) if(e[j].cap){
				if(cont[e[j].to] == 0) {
					cont[e[j].to] = 1;
					q[lq++] = e[j].to;
				}
			}
		}
		for(int i = 0; i < (n - 2) / 2; i++)
			if(cont[i] && !cont[i + (n - 2) / 2]){
				bool flag = 0;
				for(int j = head[i]; j; j = e[j].next)
					if(e[j].to == i + (n - 2) / 2)
						flag = 1;
				if(flag)
					ans.push_back(i);
			}

		return ans;
	}
} g;

const int maxn = 55;
int n;
int dis[maxn];
bool mp[maxn][maxn];
int l[maxn];
int used[maxn];

int getdis(int no) {
	int ans = 0;
	if(dis[no] != -1)
		return dis[no];
	for(int i = 0; i < n; i++) if(mp[no][i]){
		ans = max(getdis(i), ans);
	}
	return dis[no] = ans + l[no];
}

struct Farmville{
	int minTime(vector <string> s, vector <int> time, vector <int> cost, int budget){
		n = s.size();
		for(int i = 0; i < n; i++) {
			l[i] = time[i];
			for(int j = 0; j < n; j++)
				mp[i][j] = s[i][j] - '0';
		}
		ll already = 0;
		n = s.size();
		int src = n * 2, sink = n * 2 + 1;
		while(already <= budget) {
			ll maxd = 0;
			memset(dis, -1, sizeof(dis));
			for(int i = 0; i < n; i++) if(dis[i] == -1)
				getdis(i);
			memset(used, 0, sizeof(used));
			for(int i = 0; i < n; i++)
				maxd = max(maxd, 1ll * dis[i]);
			int q[maxn], lq = 0;

			g.pre(n * 2 + 2, src, sink);
			for(int i = 0; i < n; i++) {
				bool flag = 1;
				for(int j = 0; j < n; j++) if(mp[i][j])
					flag = 0;
				if(flag)
					g.addedge(src, i, INF);
			}
			for(int i = 0; i < n; i++) if(dis[i] == maxd) {
				used[i] = 1;
				q[lq++] = i;
				g.addedge(i + n, sink, INF);
			}

			for(int i = 0; i < lq; i++) {
				int no = q[i];
				if(l[no])
					g.addedge(no, no + n, cost[no]);
				else
					g.addedge(no, no + n, INF);
				for(int j = 0; j < n; j++) if(mp[no][j] && dis[j] == dis[no] - l[no]) {
					g.addedge(j + n, no, INF);
					if(used[j] == 0) {
						used[j] = 1;
						q[lq++] = j;
					}
				}
			}
			g.SAP();
			already += g.maxFlow;
			if(already > budget)
				return maxd;
			vi arr = g.getcut();
			for(int i = 0; i < (int)arr.size(); i++)
				l[arr[i]]--;
		}
		return -1;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

