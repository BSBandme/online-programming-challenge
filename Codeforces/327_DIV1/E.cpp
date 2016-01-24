//#define DEBUG
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

//待排序的字符串放在r 数组中，从r[0]到r[n-1]，长度为n，且单个字符最大值小于m。sa为结果
#ifndef DEBUG
#define maxn 10001000
#else
#define maxn  10000
#endif

const int maxl = 2;
struct trie {
	int nxt[maxl];
	int val, fail, fa;
	void init() {
		memset(nxt, -1, sizeof(nxt));
		val = -1;
		fa = -1;
	}
} arr[maxn];
int larr;
map <char, int> mp;
inline void initmp() {
	mp['a'] = 0;
	mp['b'] = 1;
}
// 返回单词最后字母的那个节点的编号，处理val需在外面处理
inline int addword(char *str, int l) {
	int no = 0;
	for(int i = 0; i < l; i++) {
		int nxt = arr[no].nxt[mp[str[i]]];
		if(nxt == -1) {
			arr[larr].init();
			arr[larr].fa = no;
			nxt = arr[no].nxt[mp[str[i]]] = larr++;
		}
		no = nxt;
	}
	return no;
}
inline void addval(int no1, int no2) {
	if(arr[no1].val == -1)
		arr[no1].val = arr[no2].val;
}

pii q[maxn];
inline void acbfs() {
	int lq = 0, rq = 0;
	for(int i = 0; i < maxl; i++) {
		if(arr[0].nxt[i] == -1) arr[0].nxt[i] = 0;
		else q[rq++] = mpr(arr[0].nxt[i], 0);
	}
	for(; lq != rq; lq++) {
		int now = q[lq].first, fail = q[lq].second;
		arr[now].fail = fail;
		for(int i = 0; i < maxl; i++) {
			if(arr[now].nxt[i] == -1) arr[now].nxt[i] = arr[fail].nxt[i];
			else {
				q[rq++] = mpr(arr[now].nxt[i], arr[fail].nxt[i]);
				addval(arr[now].nxt[i], arr[fail].nxt[i]);
			}
		}
	}
}
char orig[maxn];
int n, l;

const int maxm = 800;
int rl[maxm], last[maxm];
pii rbe[maxm];
int mat[maxm][maxm], match[maxn], bematch[maxn];
int have[maxm];

int dfs(int no) {
	if(have[no]) return 0;
	have[no] = 1;
	for(int j = 0; j < maxm; j++) if(mat[no][j] && j != bematch[no]) {
		if(bematch[j] == -1)  {
			bematch[j] = no;
			match[no] = j;
			return 1;
		} else {
			if(dfs(bematch[j])) {
				bematch[j] = no;
				match[no] = j;
				return 1;
			}
		}
	}
	return 0;
}

int ans[maxm], lans;


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
	larr = 1;
	arr[0].init();
	initmp();
	for(int i = 0; i < n; i++) {
		scanf("%s", orig);
		rl[i] = strlen(orig);
		int loc = addword(orig, rl[i]);
		arr[loc].val = i;
		last[i] = loc;
	}
	acbfs();
	for(int i = 0; i < n; i++) {
		if(arr[arr[last[i]].fail].val != -1)
			mat[arr[arr[last[i]].fail].val][i] = 1;
		for(int j = arr[last[i]].fa; j; j = arr[j].fa)
			if(arr[j].val != -1)
				mat[arr[j].val][i] = 1;
	}
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			for(int k = 0; k < n; k++)
				mat[j][k] |= mat[j][i] & mat[i][k];
	int cnt = 0;
	memset(match, -1, sizeof(match));
	memset(bematch, -1, sizeof(bematch));
	for(int i = 0; i < n; i++) {
		memset(have, 0, sizeof(have));
		if(!dfs(i)) cnt++;
	}
	for(int i = 0; i < n; i++) if(bematch[i] == -1) {
		ans[lans++] = i;
	}
	for( ; ; ) {
		bool flag = 1;
		for(int i = 0; i < lans; i++)
			for(int j = 0; j < lans; j++)
				if(mat[ans[j]][ans[i]]) {
					ans[j] = match[ans[j]];
					flag = 0;
				}
		if(flag) break;
	}
	printf("%d\n", lans);
	for(int i = 0; i < lans; i++)
		printf("%d ", ans[i] + 1);
	puts("");

    return 0;
}
