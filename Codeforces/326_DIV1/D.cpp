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

const int maxn = 50010;
pair <pii, pii> ie[maxn];
mpii du[maxn];
int n, m;
vector <pair <pii, pii> > e[maxn];
set <pii> have;
vi cir[maxn], link[maxn];
pii head[maxn];
pii val[maxn];
int nlink, ncir;
int used[maxn];
int cirans[maxn], linkans[maxn];

struct edge {
	int to, nxt;
} re[maxn * 2];
int rdu[maxn];
int rhead[maxn], le;
void addedge(int a, int b) {
	re[le].to = b;
	re[le].nxt = rhead[a];
	rhead[a] = le++;
}

void dfs(int no, int col, int fa, vi &tvi) {
	if(have.find(mpr(no, col)) != have.end())
		return;
	have.insert(mpr(no, col));
	auto p = lower_bound(e[no].begin(), e[no].end(), mpr(mpr(col, -1), mpr(-1, -1)));
 	bool flag = 0;
	for(; p != e[no].end() && p->first.first == col; p++) if(p->second.second != fa) {
		flag = 1;
		tvi.push_back(p->second.second);
		dfs(p->second.first, col, p->second.second, tvi);
		break;
	}
	if(!flag)
		head[nlink].second = no;
}

void dfs1(int no, int l, int fa) {
	bool flag = 0;
	if(l % 2 == 0)
		linkans[no] = -1;
	for(int i = rhead[no]; i != -1; i = re[i].nxt) if(re[i].nxt != fa) {
		flag = 1;
		dfs1(re[i].to, l + 1, no);
		if(l % 2 == 0)
			linkans[no] = re[i].to - nlink;
	}
}

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
#define MAXM 500000
#define MAXN 200000
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


int justdoit(int mid) {
	int rn = n + nlink;
	g.pre(rn + 2, rn, rn + 1);
	int rpn = 0;
	for(int i = 0; i < nlink; i++) {
		if(link[i].size() == 1) continue;
		if(link[i].size() & 1) {
			rpn += 2;
			if(!used[head[i].first] && !used[head[i].second] && val[i].first <= mid) {
				g.addedge(i, head[i].first + nlink, 1);
				g.addedge(i, head[i].second + nlink, 1);
			}
			if(val[i].second <= mid) {
				g.addedge(i, rn + 1, 2);
			}
			g.addedge(rn, i, 2);
		} else {
			rpn++;
			if(!used[head[i].first] && val[i].first <= mid) {
				g.addedge(i, head[i].first + nlink, 1);
			}
			if(!used[head[i].second] && val[i].second <= mid) {
				g.addedge(i, head[i].second + nlink, 1);
			}
			g.addedge(rn, i, 1);
		}
	}
	for(int i = 0; i < n; i++)
		g.addedge(i + nlink, rn + 1, 1);
	g.SAP();
	return g.maxFlow >= rpn;
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

	scanf("%d%d", &n, &m);
	for(int i = 0; i < m; i++) {
		int u, v, c, t;
		scanf("%d%d%d%d", &u, &v, &c, &t);
		u--; v--;
		e[v].push_back(mpr(mpr(c, t), mpr(u, i)));
		e[u].push_back(mpr(mpr(c, t), mpr(v, i)));
		ie[i] = mpr(mpr(c, t), mpr(u, v));
		du[u][c]++;
		du[v][c]++;
	}
	for(int i = 0; i < n; i++)
		sort(e[i].begin(), e[i].end());
	bool flag =  1;
	for(int i = 0; i < n; i++)
		for(auto it : du[i])
			if(it.second > 2)
				flag = 0;

	if(!flag) {
		puts("No");
		return 0;
	}

	for(int i = 0; i < n; i++) {
		for(auto it : du[i]) if(it.second == 1 && have.find(mpr(i, it.first)) == have.end()) {
			head[nlink].first = i;
			dfs(i, it.first, -1, link[nlink]);
			nlink++;
		}
	}
	for(int i = 0; i < n; i++) {
		for(auto it : du[i]) if(have.find(mpr(i, it.first)) == have.end()) {
			dfs(i, it.first, -1, cir[ncir]);
			ncir++;
		}
	}
	for(int i = 0; i < nlink; i++) {
		for(auto it : link[i]) {
			int u = ie[it].second.first, v = ie[it].second.second;
			if(u != head[i].first && u != head[i].second)
				used[u]++;
			if(v != head[i].first && v != head[i].second)
				used[v]++;
		}
	}
	int ans = 0;
	for(int i = 0; i < ncir; i++) {
		for(auto it : cir[i]) {
			int u = ie[it].second.first, v = ie[it].second.second;
			used[u]++;
			used[v]++;
		}
		if(cir[i].size() & 1)
			flag = 0;
		int rans = 0, rans1 = 0;
		for(int j = 0; j < (int)cir[i].size(); j += 2)
			rans = max(rans, ie[cir[i][j]].first.second);
		for(int j = 1; j < (int)cir[i].size(); j += 2)
			rans1 = max(rans1, ie[cir[i][j]].first.second);
		if(rans > rans1)
			cirans[i] = 1;
		else
			cirans[i] = 0;
		ans = max(ans, min(rans, rans1));
	}
	for(int i = 0; i < n; i++)
		if(used[i] > 2 || (used[i] & 1))
			flag = 0;
	if(!flag) {
		puts("No");
		return 0;
	}

	for(int i = 0; i < nlink; i++) {
		int rans = 0;
		for(int j = 0; j < (int)link[i].size(); j += 2) {
			int eno = link[i][j];
			rans = max(ie[eno].first.second, rans);
		}
		val[i].first = rans;
		rans = 0;
		for(int j = 1; j < (int)link[i].size(); j += 2) {
			int eno = link[i][j];
			rans = max(ie[eno].first.second, rans);
		}
		val[i].second = rans;
	}

//	if(n == 33458) {
//		cout << ans << ' ' << ncir << ' ' << nlink << endl;
//		int p = 0;
//		int rmax = 0;
//		for(int i = 0; i < nlink; i++) {
//			if(link[i].size() > 1) {
//				p++;
//				rmax = max(rmax, min(val[i].first, val[i].second));
//				cout << head[i].first << ' ' << head[i].second << ' ';
//			}
//		}
//		cout << p << ' ' << rmax << endl;
//	}

	int be = 0, en = MOD;
	int rn = nlink + n;
	for(; be != en && flag; ) {
		int mid = (be + en) / 2;
		if(justdoit(mid)) en = mid;
		else be = mid + 1;
	}

	if(!justdoit(be)) flag = 0;
	if(!flag) {
		puts("No");
		return 0;
	}
	for(int i = 0; i < nlink; i++)
		linkans[i] = -1;
	for(int i = 2; i < g.p; i += 2) {
		if(g.e[i + 1].to < nlink) {
			int no = g.e[i + 1].to;
			int rno = g.e[i].to - nlink;
			if(g.e[i].to == rn) continue;
			if(g.e[i].to == rn + 1 && g.e[i].cap != 2)
				linkans[no] = rno;
			else if(g.e[i].cap == 0) {
				linkans[no] = rno;
			}
		}
	}
	puts("Yes");
	vi anse;
	for(int i = 0; i < nlink; i++) {
		if(link[i].size() == 1)
			continue;
		if(linkans[i] == head[i].first || (link[i].size() % 2 && linkans[i] == head[i].second)) {
			ans = max(ans, val[i].first);
			for(int j = 0; j < (int)link[i].size(); j += 2)
				anse.push_back(link[i][j]);
		} else {
			ans = max(ans, val[i].second);
			for(int j = 1; j < (int)link[i].size(); j += 2)
				anse.push_back(link[i][j]);
		}
	}
	for(int i = 0; i < ncir; i++) {
		for(int j = cirans[i]; j < (int)cir[i].size(); j += 2)
			anse.push_back(cir[i][j]);
	}
	printf("%d %d\n", ans, (int)anse.size());
	for(int i = 0; i < (int)anse.size(); i++)
		printf("%d ", anse[i] + 1);

//	for(int i = 0; i < nlink; i++)
//		linkans[i] = -2;
//	int q[maxn], lq = 0;
//	for(int i = 0; i < nlink; i++) if(rdu[i] == 1)
//		q[lq++] = i;
//	for(int i = 0; i < lq; i++) {
//		int no = q[i];
//		int rno = -1;
//		for(int j = rhead[no]; j != -1; j = re[j].nxt) if(!used[re[j].to]) {
//			rno = re[j].to;
//			used[re[j].to] = 1;
//		}
//		linkans[no] = rno;
//		for(int j = rhead[rno]; j != -1; j = re[j].nxt) {
//			rdu[re[j].to]--;
//			if(rdu[re[j].to] == 1)
//				q[lq++] = re[j].to;
//		}
//	}
//	for(int i = 0; i < nlink; i++) if(linkans[i] == -2) {
//
//	}



    return 0;
}
