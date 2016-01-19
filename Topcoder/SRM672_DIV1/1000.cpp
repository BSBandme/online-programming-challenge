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

const int maxn = 55;
const int maxh = 11;

ll canreach[maxn][maxh];
int mp[maxn][maxn];

struct node {
	ll have, no, rh, dis;
	node(ll a = 0, ll b = 0, ll c = 0, ll d = 0) :
		have(a), no(b), rh(c), dis(d) {}
	bool operator < (const node &a) const {
		return dis < a.dis;
	}
	bool operator > (const node &a) const {
		return dis > a.dis;
	}
};
priority_queue <node, vector <node>, greater <node> > pq;
int n;
int ansh[maxn][maxh];
int tmp[maxn][maxh][maxn];

struct Tdetective{
	int reveal(vector <string> s){
		n = s.size();
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++)
				mp[i][j] = s[i][j] - '0';
		}

		memset(canreach, 0, sizeof(canreach));
		memset(ansh, 0x3f, sizeof(ansh));
		for(int h = 9; h >= 0; h--) {
			for(int i = 0; i < n; i++)
				canreach[i][h] |= 1ll << i;
			int rmp[maxn];
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++)
					rmp[j] = mp[i][j];
				ll searched = 0;
				for(bool flag = 1; flag;) {
					flag = 0;
					for(int j = 0; j < n; j++) if(rmp[j] >= h && (searched & (1ll << j)) == 0) {
						flag = 1;
						ll rmask = canreach[j][rmp[j]] ^ (canreach[j][rmp[j]] & canreach[i][h]);
						for(int p = 0; p < n; p++) if(rmask & (1ll << p))
							for(int k = 0; k < n; k++)
								rmp[k] = max(rmp[k], mp[p][k]);
						canreach[i][h] |= canreach[j][rmp[j]];
						searched |= 1ll << j;
					}
				}
			}
		}

		memset(tmp, 0, sizeof(tmp));
		pq.push(node(1, 0, 0, 0));
		for(int j = 0; j < maxh; j++)
			ansh[0][j] = 0;
		for(int j = 0; j < n; j++) {
			tmp[0][0][j] = mp[0][j];
		}
		for(; pq.size(); ) {
			ll rhave = pq.top().have;
			ll no = pq.top().no;
			ll dis = pq.top().dis;
			ll h = pq.top().rh;
			pq.pop();
			if(ansh[no][h] != dis)
				continue;
			for(int j = 9; j >= h; j--) {
				ll phave = rhave;
				int rmp[maxn];
				memcpy(rmp, tmp[no][h], sizeof(rmp));
				for(int i = 0; i < n; i++)
					if(tmp[no][h][i] == j && (rhave & (1ll << i)) == 0) {
 						ll rmask = canreach[i][tmp[no][h][i]] ^ (canreach[i][tmp[no][h][i]] & phave);
						for(int j = 0; j < n; j++) if(rmask & (1ll << j))
							for(int k = 0; k < n; k++)
								rmp[k] = max(rmp[k], mp[j][k]);
						phave |= canreach[i][tmp[no][h][i]];
						if(ansh[i][tmp[no][h][i]] > bitnum(rhave)) {
							ansh[i][tmp[no][h][i]] = bitnum(rhave);
							pq.push(node(rhave | 1ll << i, i, tmp[no][h][i], bitnum(rhave)));
							for(int ii = 0; ii < n; ii++) if((rhave & (1ll << ii)) || ii == i)
								for(int jj = 0; jj < n; jj++) {
									tmp[i][tmp[no][h][i]][jj] = max(tmp[i][tmp[no][h][i]][jj], mp[ii][jj]);
								}
						}
					}
				memcpy(tmp[no][h], rmp, sizeof(rmp));
				rhave = phave;
			}
		}
		ll pans = 0;
		for(int i = 1; i < n; i++) {
			int ans = MOD, ansno = -1;
			for(int j = 0; j < 10; j++) if(ansh[i][j] < ans) {
				ans = ansh[i][j];
				ansno = j;
			}
//			cerr << i << "|" << ansno << "|" << ans << endl;
			pans += i * ans;
		}

		return pans;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

