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

const int maxn = 1300;
const int maxp = 8;
int dp1[2][maxn][maxp][maxp];
int dp2[2][maxn][maxp][maxp];
int n;

struct BearPairs{
	int minCost(string str, vector <int> cost, int rk){
		n = str.size();
		memset(dp1, 0x3f, sizeof(dp1));
		memset(dp2, 0x3f, sizeof(dp2));
		for(int i = 0; i < 6; i++) {
			dp2[0][0][i][0] = 0;
			dp1[0][0][i][0] = 0;
		}

		int now = 0, nxt = 1;
		for(int i = 0; i < n; i++) {
			memset(dp1[nxt], 0x3f, sizeof(dp1[nxt]));
			memset(dp2[nxt], 0x3f, sizeof(dp2[nxt]));

			// important!!!------------------------------------
			// "aaaaccbbbbbabaacac", rk == 0, costs all == 0, the answer is 2700
			for(int j = 0; j <= min(i, n - i); j++) for(int k = 0; k < 6; k++) for(int t = 0; t <= rk; t++) if(dp1[now][j][k][t] < MOD) {
				int val = dp1[now][j][k][t];
				for(int ri = 0; ri < 6; ri++) if(ri != k)
					dp2[now][j][ri][t] = min(dp2[now][j][ri][t], val);
			}
			// important!!!------------------------------------


			for(int j = 0; j <= min(i, n - i); j++) for(int k = 0; k < 6; k++) for(int t = 0; t <= rk; t++) if(dp1[now][j][k][t] < MOD) {
				dp1[nxt][j][k][t + 1] = min(dp1[now][j][k][t], dp1[nxt][j][k][t + 1]);
				if(j == 0) {
					for(int ri = 0; ri < 6; ri++) if(ri != str[i] - 'a')
						dp2[nxt][1][ri][t] = min(dp1[now][j][k][t] + cost[i] - 100 * i, dp2[nxt][1][ri][t]);
					dp1[nxt][1][str[i] - 'a'][t] = min(dp1[now][j][k][t] + cost[i] - 100 * i, dp1[nxt][1][str[i] - 'a'][t]);
				} else if(str[i] - 'a' == k){
					dp1[nxt][j + 1][k][t] = min(dp1[now][j][k][t] + cost[i] - 100 * i, dp1[nxt][j + 1][k][t]);
				} else {
					dp1[nxt][j - 1][k][t] = min(dp1[now][j][k][t] + cost[i] + 100 * i, dp1[nxt][j - 1][k][t]);
				}
			}
			for(int j = 0; j <= min(i, n - i); j++) for(int k = 0; k < 6; k++) for(int t = 0; t <= rk; t++) if(dp2[now][j][k][t] < MOD) {
				int val = dp2[now][j][k][t];
				dp2[nxt][j][k][t + 1] = min(val, dp2[nxt][j][k][t + 1]);
				if(j == 0) {
					for(int ri = 0; ri < 6; ri++) if(ri != str[i] - 'a')
						dp2[nxt][1][ri][t] = min(val + cost[i] - 100 * i, dp2[nxt][1][ri][t]);
					dp1[nxt][1][str[i] - 'a'][t] = min(val + cost[i] - 100 * i, dp1[nxt][1][str[i] - 'a'][t]);
				} else if(str[i] - 'a' == k) {
					int &rval = dp2[nxt][j - 1][k][t];
					rval = min(rval, val + cost[i] + 100 * i);
				} else {
					int &rval = dp2[nxt][j + 1][k][t];
					rval = min(rval, val + cost[i] - 100 * i);
				}
			}
			swap(now, nxt);
		}

		int ans = MOD * 2;
		for(int i = 0; i < 6; i++) for(int j = 0; j <= rk; j++)
			ans = min1(dp1[now][0][i][j], dp2[now][0][i][j], ans);
		if(ans > MOD) ans = -1;
		return ans;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); if ((Case == -1) || (Case == 5)) test_case_5(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const int &Expected, const int &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() {
		string Arg0 = "aaaaccbbbbbabaacac";
		int Arr1[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		vector<int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0])));
		int Arg2 = 0;
		int Arg3 = 20;
		verify_case(0, Arg3, minCost(Arg0, Arg1, Arg2));
	}
	void test_case_1() { string Arg0 = "cdbcadc"; int Arr1[] = {261,208,150,250,92,226,176}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arg2 = 1; int Arg3 = 1402; verify_case(1, Arg3, minCost(Arg0, Arg1, Arg2)); }
	void test_case_2() { string Arg0 = "deebaffafdaaceaa"; int Arr1[] = {160,268,253,210,34,28,180,70,5,42,177,234,108,117,215,1}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arg2 = 2; int Arg3 = 2507; verify_case(2, Arg3, minCost(Arg0, Arg1, Arg2)); }
	void test_case_3() { string Arg0 = "babbbabbbbababababbb"; int Arr1[] = {184,189,202,170,296,71,136,48,51,161,221,24,221,186,223,228,73,274,279,22}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arg2 = 4; int Arg3 = -1; verify_case(3, Arg3, minCost(Arg0, Arg1, Arg2)); }
	void test_case_4() { string Arg0 = "aaaaaaaaaaaaaaaaaa"; int Arr1[] = {237,185,24,175,107,251,299,81,282,20,150,164,240,225,166,261,164,123}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arg2 = 4; int Arg3 = -1; verify_case(4, Arg3, minCost(Arg0, Arg1, Arg2)); }
	void test_case_5() { string Arg0 = "acadeffbffbfccbe"; int Arr1[] = {62,113,189,161,211,150,69,60,99,212,37,274,110,265,199,192}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arg2 = 4; int Arg3 = 2235; verify_case(5, Arg3, minCost(Arg0, Arg1, Arg2)); }

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	BearPairs ___test;
	___test.run_test(-1);

	return 0;
}
// END CUT HERE


