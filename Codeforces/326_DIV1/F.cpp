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

//待排序的字符串放在r 数组中，从r[0]到r[n-1]，长度为n，且单个字符最大值小于m。sa为结果
#define maxn 200100
int wa[maxn * 2],wb[maxn * 2],wv[maxn * 2],ws1[maxn * 2];
int cmp(int *r,int a,int b,int l)
{return r[a]==r[b]&&r[a+l]==r[b+l];}
void da(int *r,int *sa,int n,int m) {
	 int i,j,p,*x=wa,*y=wb,*t;
	 for(i=0;i<m;i++) ws1[i]=0;
	 for(i=0;i<n;i++) ws1[x[i]=r[i]]++;
	 for(i=1;i<m;i++) ws1[i]+=ws1[i-1];
	 for(i=n-1;i>=0;i--) sa[--ws1[x[i]]]=i;
	 for(j=1,p=1;p<=n;j*=2,m=p) {
	   for(p=0,i=n-j;i<n;i++) y[p++]=i;
	   for(i=0;i<n;i++) if(sa[i]>=j) y[p++]=sa[i]-j;
	   for(i=0;i<n;i++) wv[i]=x[y[i]];
	   for(i=0;i<m;i++) ws1[i]=0;
	   for(i=0;i<n;i++) ws1[wv[i]]++;
	   for(i=1;i<m;i++) ws1[i]+=ws1[i-1];
	   for(i=n-1;i>=0;i--) sa[--ws1[wv[i]]]=y[i];
	   for(t=x,x=y,y=t,p=2,x[sa[0]]=1,i=1;i<n;i++)
	   x[sa[i]]=cmp(y,sa[i-1],sa[i],j)?p-1:p++;
	 }
	 return;
}
void calheight(int *r, int *sa, int*rank, int*height,int n) {
	 int i,j,k=0;
	 for(i=0; i<n;i++) rank[sa[i]]=i;
	 for(i=0;i<n;i++){
		 if(rank[i] == 0) continue;
		 for(k?k--:0,j=sa[rank[i]-1]; i + k < n && j + k < n && r[i+k]==r[j+k];k++);
		 height[rank[i] - 1]=k;
	 }
	 return;
}
int d[20];
int logg[maxn];
int st[maxn][22];
// h为height数组，n为长度
void InitRMQ(const int &n, int *h) {
	int i, j;
	for( d[0]=1, i=1; i < 21; ++i ) d[i] = 2*d[i-1];
	for( i=0; i < n; ++i ) st[i][0] = h[i];
	int k = int( log(double(n))/log(2.0) ) + 1;
	logg[0] = 0; logg[1] = 0;
	for(int i = 2; i <= n; i++) logg[i] = logg[i / 2] + 1;
	for( j=1; j < k; ++j )
		for( i=0; i < n; ++i ) {
			if( i+d[j-1]-1 < n ){
				st[i][j] = min(st[i][j-1],
						st[i+d[j-1]][j-1]);
			}
			else break; // st[i][j] = st[i][j-1];
		}
}
//x, y为rank的序号(在sa中为x和y位置的lcp应为quelcp(x, y))
inline int quelcp(int x, int y) {
	int k = logg[y - x + 1];
	return min(st[x][k], st[y-d[k]+1][k]);
}
//此为真正的quelcp(y--）
inline int quelcp1(int x, int y) {
	if(y < x) swap(x, y);
	y--;
	int k = logg[y - x + 1];
	return min(st[x][k], st[y-d[k]+1][k]);
}
//querange为询问rank序号为no的位置与之lcp >= height的左右边界
inline pii querange(int no, int height, int l) {
	int be1 = 0, en1 = no - 1;
	while(be1 < en1) {
		int height1 = (be1 + en1) / 2;
		if(quelcp(height1, no - 1) < height) be1 = height1 + 1;
		else en1 = height1;
	}
	if(quelcp(be1, no - 1) < height) be1++;
	int from = min1(be1, no);
	be1 = no; en1 = l - 1;
	while(be1 < en1) {
		int height1 = (be1 + en1 + 1) / 2;
		if(quelcp(no, height1) < height) en1 = height1 - 1;
		else be1 = height1;
	}
	if(quelcp(no, be1) < height) be1--;
	int to = max1(be1 + 1, no);
	return mpr(from, to);
}

int arr[maxn], sa[maxn], rk[maxn], h[maxn];
int n, nq;
int larr;
char tstr[maxn];
int be[maxn], len[maxn], belong[maxn];
pii range[maxn];
int cnt[maxn];
ll all[maxn];
ll ans[maxn];
int que[maxn][3];
int yu;
vi stat[maxn / 2];
pii already[maxn];
vi hq[maxn];

int main() {

//............................不要再忘了检查maxn大小了！！！！BSBandme你个SB！！！！...................................................

	ios_base::sync_with_stdio(0);
	#ifdef DEBUG //......................................................................................................
	freopen("input.txt", "r", stdin);
	int __size__ = 256 << 20; // 256MB
	char *__p__ = (char*)malloc(__size__) + __size__;
	__asm__("movl %0, %%esp\n" :: "r"(__p__));
	#endif //...........................................................................................................

	scanf("%d%d", &n, &nq);
	for(int i = 0; i < n; i++) {
		scanf("%s", tstr);
		be[i] = larr;
		int ltstr = strlen(tstr);
		len[i] = ltstr;
		for(int j = 0; j < ltstr; j++) {
			belong[larr] = i;
			arr[larr++] = tstr[j] - 'a' + 1;
		}
		belong[larr] = -1;
		arr[larr++] = i + 26 + 1;
	}
	da(arr, sa, larr, n + 30 + 26);
	calheight(arr, sa, rk, h, larr);
	InitRMQ(larr, h);
	for(int i = 0; i < n; i++)
		range[i] = querange(rk[be[i]], len[i], larr);

	for(int i = 0; i < nq; i++) {
		for(int j = 0; j < 3; j++) {
			scanf("%d", que[i] + j);
			que[i][j]--;
		}
		hq[que[i][2]].emplace_back(i);
		already[i] = mpr(larr, -1);
	}

	yu = 150;
	for(int i = 0; i < n; i += yu) {
		int en = min(n - 1, i + yu - 1);
		for(int i = 0; i < larr; i++) {
			all[i] = 0;
			cnt[i] = 0;
		}
		for(int j = i; j <= en; j++) {
			cnt[range[j].first]++;
			cnt[range[j].second + 1]--;
		}
		for(int j = 0, rcnt = 0; j < larr; j++) {
			rcnt += cnt[j];
			if(belong[sa[j]] == -1) break;
			all[belong[sa[j]]] += rcnt;
		}
		for(int j = 0; j < nq; j++)
			if(que[j][0] <= i && que[j][1] >= en) {
				ans[j] += all[que[j][2]];
				already[j].first = min(already[j].first, i);
				already[j].second = max(already[j].second, en);
			}
	}
	memset(cnt, 0, sizeof(cnt));
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < (int)hq[i].size(); j++) {
			int qno = hq[i][j];
			if(already[qno].first > already[qno].second) {
				for(int k = que[qno][0]; k <= que[qno][1]; k++)
					stat[range[k].first].emplace_back(qno);
			} else {
				for(int k = que[qno][0]; k < already[qno].first; k++)
					stat[range[k].first].emplace_back(qno);
				for(int k = already[qno].second + 1; k <= que[qno][1]; k++)
					stat[range[k].first].emplace_back(qno);
			}
		}
	}
	for(int i = 0; i < larr; i++) {
		for(int j = 0; j < (int)stat[i].size(); j++) {
			int qno = stat[i][j];
			ans[qno] -= cnt[que[qno][2]];
		}
		cnt[belong[sa[i]]]++;
		if(belong[sa[i]] == -1) break;
		stat[i].clear();
	}
	memset(cnt, 0, sizeof(cnt));
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < (int)hq[i].size(); j++) {
			int qno = hq[i][j];
			if(already[qno].first > already[qno].second) {
				for(int k = que[qno][0]; k <= que[qno][1]; k++)
					stat[range[k].second + 1].emplace_back(qno);
			} else {
				for(int k = que[qno][0]; k < already[qno].first; k++)
					stat[range[k].second + 1].emplace_back(qno);
				for(int k = already[qno].second + 1; k <= que[qno][1]; k++)
					stat[range[k].second + 1].emplace_back(qno);
			}
		}
	}
	for(int i = 0; i < larr; i++) {
		for(int j = 0; j < (int)stat[i].size(); j++) {
			int qno = stat[i][j];
			ans[qno] += cnt[que[qno][2]];
		}
		cnt[belong[sa[i]]]++;
		if(belong[sa[i]] == -1) break;
		stat[i].clear();
	}
	for(int i = 0; i < nq; i++)
		printf("%I64d\n", ans[i]);


    return 0;
}
