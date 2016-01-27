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

const int mod = MOD;
const int maxn = 55;
const int maxn1 = 50000;
const int maxl = 10;
int d;
string arc, be, en;
int dp[maxn][2][maxn1];

struct trie {
	int nxt[maxl];
	int deep;
	int val;
	void init() {
		memset(nxt, -1, sizeof(nxt));
		val = 0;
		deep = 0;
	}
} arr[maxn1];
int larr;
map <char, int> mp;
inline void initmp() {
	for(char ch = '0'; ch <= '9'; ch++)
		mp[ch] = ch - '0';
}
// 返回单词最后字母的那个节点的编号，处理val需在外面处理
inline int addword(string str) {
	int l = str.size(), no = 0;
	for(int i = 0; i < l; i++) {
		int nxt = arr[no].nxt[mp[str[i]]];
		if(nxt == -1) {
			arr[larr].init();
			nxt = arr[no].nxt[mp[str[i]]] = larr++;
		}
		arr[nxt].deep = i + 1;
		no = nxt;
	}
	return no;
}
inline void addval(int no1, int no2) {
	arr[no1].val |= arr[no2].val;
}
inline void acbfs() {
	pii q[maxn1];
	int lq = 0, rq = 0;
	for(int i = 0; i < maxl; i++) {
		if(arr[0].nxt[i] == -1) arr[0].nxt[i] = 0;
		else q[rq++] = mpr(arr[0].nxt[i], 0);
	}
	for(; lq != rq; lq++) {
		int now = q[lq].first, fail = q[lq].second;
		for(int i = 0; i < maxl; i++) {
			if(arr[now].nxt[i] == -1) arr[now].nxt[i] = arr[fail].nxt[i];
			else {
				q[rq++] = mpr(arr[now].nxt[i], arr[fail].nxt[i]);
				addval(arr[now].nxt[i], arr[fail].nxt[i]);
			}
		}
	}
}


ll getdp(string str, int rl, int cate) {
	for(int i = 0; i < (int)str.size() + 3; i++)
		for(int j = 0; j < 2; j++)
			for(int k = 0; k < larr; k++) {
				dp[i][j][k] = 0;
			}
	dp[0][0][0] = 1;
	for(int i = 0; i < (int)str.size(); i++) {
		for(int j = 0; j < 2; j++) {
			int maxv = str[i] - '0';
			if(j) maxv = 9;
			for(int k = 0; k < larr; k++) if(dp[i][j][k]) {
				for(int t = 0; t <= maxv; t++) {
					int p = arr[k].nxt[t];
					if(arr[k].deep == rl)
						p = k;
					add(dp[i + 1][j | (t != maxv)][p], dp[i][j][k]);
				}
			}
		}
	}
	ll ans = 0;
	for(int i = 0; i < larr; i++) if(arr[i].deep == rl)
		ans += cate * dp[(int)str.size()][0][i] + dp[(int)str.size()][1][i];
	return ans % mod;
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

	cin >> arc >> be >> en;
	int larc = arc.size();
	d = be.size() / 2;
	larr = 1;
	arr[0].init();
	initmp();
	for(int i = 0; i <= larc - d; i++)
		addword(arc.substr(i, d));
	acbfs();

//	for(int i = d; i <= (int)be.size(); i++) {
//		for(int j = 0; j <= larc - i; j++) {
//			string tstr = arc.substr(j, i);
//			if(have.find(tstr) != have.end())
//				continue;
//			kmppre(tstr);
//			ll tans = getdp(en, i, tstr, 1) - getdp(be, i, tstr, 0);
//			have[tstr].first = tans % mod;
//			set <string> rs;
//			for(int k = d; k < i; k++) {
//				for(int t = 0; t <= i - k; t++)
//					rs.insert(tstr.substr(t, k));
//			}
//			for(auto it : rs) {
//				have[tstr].second += have[it].second;
//			}
//			have[tstr].second = 1 - have[tstr].second;
//			ans += have[tstr].second * have[tstr].first % mod;
//		}
//	}
	cout << (mod + getdp(en, d, 1) - getdp(be, d, 0)) % mod << endl;

    return 0;
}
