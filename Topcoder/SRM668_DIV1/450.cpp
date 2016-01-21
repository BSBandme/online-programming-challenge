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

const int maxn = 4010;
struct edge {
	int to, nxt;
} e[maxn * 2];
int head[maxn], le;
int n, m;
int dp[maxn][maxn];
int q[maxn * maxn][2], lq;

bool can(int from, int to) {
	memset(dp, 0, sizeof(dp));
	dp[from][0] = 1;
	lq = 0;
	q[lq][0] = from;
	q[lq++][1] = 0;
	for(int i = 0; i < lq; i++) {
		int no = q[i][0], cate = q[i][1];
		if(cate >= m * 2) continue;
		for(int j = head[no]; j != -1; j = e[j].nxt) {
			if(dp[e[j].to][cate + 1] == 0) {
				dp[e[j].to][cate + 1] = 1;
				q[lq][0] = e[j].to;
				q[lq++][1] = 1 + cate;
			}
		}
	}
	int ans = 0;
	bool have = 0;
	for(int j = 0; j <= m * 2; j++) {
		if(dp[from][j]) ans = __gcd(ans, j);
		if(dp[to][j]) have = 1;
	}
	return ans == 1 && have;
}

struct WalkingToSchool{
	string canWalkExactly(int N, vector <int> from, vector <int> to){
		n = N;
		m = from.size();
		memset(head, -1, sizeof(head)	);
		for(int i = 0; i < (int)from.size(); i++) {
			e[le].to = to[i];
			e[le].nxt = head[from[i]];
			head[from[i]] = le++;
		}

		if(can(0, 1) && can(1, 0))
			return "Freedom";

		return "Chores";
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); if ((Case == -1) || (Case == 5)) test_case_5(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const string &Expected, const string &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() {
		int Arg0 = 188;
		int Arr1[] = {160, 164, 184, 156, 37, 172, 166, 58, 84, 170, 16, 96, 72, 48, 150, 187, 176, 31, 160, 47, 56, 125, 79, 24, 166, 178, 169, 132, 66, 10, 36, 19, 113, 12, 154, 145, 145, 16, 110, 119, 107, 134, 141, 153, 99, 60, 138, 49, 84, 180, 78, 86, 69, 175, 44, 137, 148, 39, 170, 44, 157, 49, 97, 41, 172, 90, 132, 42, 86, 110, 13, 80, 140, 80, 162, 187, 148, 187, 108, 127, 165, 107, 115, 114, 12, 108, 11, 99, 158, 165, 151, 98, 130, 1, 167, 51, 51, 109, 44, 36, 26, 26, 145, 57, 148, 45, 9, 1, 179, 78, 97, 34, 97, 84, 157, 49, 127, 187, 76, 123, 123, 153, 71, 160, 6, 0, 134, 115, 181, 44, 14, 104, 44, 119, 187, 128, 62, 12, 82, 124, 70, 171, 152, 65, 16, 3, 153, 5, 92, 6, 94, 164, 174, 31, 63, 176, 11, 82, 68, 165, 122, 163, 140, 52, 52, 135, 72, 84, 18, 66, 176, 128, 111, 123, 47, 22, 34, 50, 137, 82, 107, 76, 27, 67, 55, 72, 131, 72, 183, 156, 99, 156, 161, 30, 166, 118, 62, 123, 130, 121, 12, 86, 186, 93, 177, 42, 25, 79, 30, 75, 73, 9, 73, 77, 142, 170, 5, 19, 41, 121, 176, 114, 159, 166, 63, 92, 61, 175, 125, 0, 63, 173, 90, 148, 187, 178, 24, 28, 53, 32, 96, 93, 31, 104, 32, 43, 47, 152, 140, 79, 112, 113, 136, 181, 48, 85, 163, 166, 157, 181, 43, 130, 10, 183, 142, 147, 98, 169, 153, 95, 53, 164, 76, 94, 4, 22, 0, 101, 141, 75, 3, 47, 0, 51, 11, 180, 147, 150, 159, 156, 174, 34, 151, 160, 95, 152, 176, 113, 98, 171, 73, 62, 158, 134, 151, 150, 125, 134, 50, 123, 116, 29, 49, 124, 4, 171, 81, 27, 70, 89, 9, 19, 84, 73, 11, 118, 60, 95, 41, 83, 159, 86, 161, 178, 175, 183, 45, 90, 84, 139, 103, 53, 92, 57, 114, 48, 34, 182, 50, 186, 2, 158, 43, 153, 139, 186, 38, 77, 187, 66, 94, 2, 59, 164, 37, 150, 185, 174, 3, 128, 3, 143, 40, 177, 179, 30, 1, 150, 62, 26, 38, 156, 2, 56, 164, 139, 45, 12, 59, 78, 75, 149, 69, 129, 179, 5, 133, 42, 166, 131, 92, 41, 17, 159, 98, 67, 35, 49, 62, 6, 41, 177, 57, 42, 20, 52, 63, 49, 2, 93, 30, 39, 36, 168, 32, 10, 57, 178, 123, 56, 138, 97, 183, 2, 166, 85, 135, 25, 83, 147, 159, 34, 28, 41, 146, 131, 55, 34, 159, 17, 157, 132, 57, 51, 66, 107, 92, 164, 146, 184, 102, 153, 50, 46, 8, 155, 36, 47, 58, 17, 160, 50, 89, 145, 91, 61, 17, 1, 143, 170, 69, 143, 28, 113, 167, 18, 95, 106, 129, 16, 177, 44, 147, 27, 106, 164, 180, 7, 69, 12, 13, 157, 81, 146, 106, 145, 85, 57, 186, 136, 12, 58, 51, 146, 45, 32, 163, 41, 131, 30, 148, 164, 174, 133, 184, 92, 168, 51, 52, 170, 48, 17, 121, 134, 178, 129, 182, 120, 107, 51, 27, 14, 16, 149, 161, 50, 6, 186, 16, 97, 177, 157, 123, 16, 180, 42, 110, 144, 172, 45, 92, 103, 144, 80, 12, 174, 25, 185, 92, 156, 98, 31, 117, 163, 125, 81, 169, 78, 13, 152, 67, 5, 94, 104, 126, 62, 123, 39, 106, 123, 38, 56, 181, 164, 90, 109, 79, 49, 151, 69, 69, 78, 17, 38, 169, 65, 97, 124, 60, 12, 174, 116, 158, 179, 153, 21, 54, 17, 5, 45, 43, 100, 95, 40, 47, 27, 139, 103, 123, 159, 130, 148, 183, 139, 18, 107, 73, 184, 99, 81, 80, 147, 19, 95, 48, 144, 45, 155, 71, 35, 83, 73, 187, 100, 19, 130, 20, 75, 157, 16, 150, 169, 86, 163, 29, 128, 2, 177, 4, 6, 159, 47, 38, 3, 131, 1, 43, 50, 154, 165, 166, 155, 104, 148, 140, 100, 75, 151, 184, 182, 175, 2, 152, 144, 179, 136, 118, 124, 142, 129, 96};
		vector<int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0])));
		int Arr2[] = {131, 70, 116, 154, 133, 50, 149, 115, 30, 81, 86, 75, 140, 75, 47, 8, 71, 85, 58, 133, 130, 110, 63, 59, 155, 168, 7, 24, 165, 66, 6, 177, 27, 96, 23, 83, 180, 108, 44, 131, 31, 33, 64, 31, 132, 67, 106, 3, 134, 76, 85, 15, 11, 110, 28, 37, 79, 187, 129, 74, 119, 163, 168, 87, 17, 115, 150, 25, 134, 157, 28, 4, 94, 156, 142, 185, 88, 162, 162, 61, 172, 114, 21, 174, 113, 38, 173, 122, 170, 108, 76, 171, 120, 130, 148, 159, 176, 46, 157, 114, 114, 4, 182, 44, 139, 22, 180, 122, 53, 83, 144, 165, 171, 9, 103, 184, 132, 181, 95, 167, 73, 127, 69, 151, 25, 9, 21, 78, 114, 111, 179, 157, 99, 69, 138, 0, 97, 149, 93, 102, 112, 0, 120, 68, 60, 130, 178, 172, 22, 168, 162, 27, 39, 172, 86, 49, 125, 16, 29, 105, 123, 148, 74, 63, 131, 103, 28, 16, 187, 146, 48, 74, 134, 154, 158, 110, 93, 48, 62, 167, 133, 128, 158, 110, 93, 148, 62, 145, 9, 16, 182, 166, 48, 14, 63, 183, 83, 21, 140, 9, 28, 131, 0, 99, 96, 28, 162, 74, 46, 130, 8, 178, 29, 96, 79, 146, 46, 143, 140, 85, 157, 172, 174, 24, 108, 70, 154, 102, 112, 149, 98, 122, 169, 14, 32, 53, 83, 40, 49, 134, 13, 71, 168, 20, 86, 96, 131, 107, 184, 138, 158, 32, 26, 161, 112, 44, 51, 178, 152, 26, 15, 151, 45, 31, 183, 174, 163, 118, 147, 155, 43, 174, 120, 76, 183, 16, 59, 2, 73, 120, 78, 14, 91, 155, 76, 2, 105, 79, 87, 92, 180, 48, 91, 132, 143, 47, 78, 79, 81, 163, 79, 16, 7, 105, 68, 159, 161, 69, 98, 5, 146, 100, 173, 21, 50, 66, 107, 130, 173, 15, 8, 148, 86, 147, 16, 62, 143, 119, 73, 102, 185, 11, 140, 66, 10, 84, 37, 59, 41, 33, 125, 134, 113, 104, 13, 148, 7, 180, 14, 3, 151, 33, 105, 29, 160, 25, 21, 149, 62, 181, 123, 20, 139, 69, 42, 60, 159, 142, 169, 99, 17, 13, 148, 176, 128, 118, 126, 0, 101, 126, 178, 108, 48, 34, 43, 112, 10, 31, 38, 75, 6, 74, 21, 126, 135, 183, 71, 116, 97, 4, 112, 177, 157, 137, 64, 137, 23, 79, 52, 155, 11, 82, 14, 136, 146, 84, 124, 106, 4, 2, 109, 172, 113, 124, 117, 86, 38, 138, 43, 115, 143, 78, 117, 168, 90, 74, 88, 177, 173, 25, 64, 148, 67, 60, 158, 12, 99, 118, 148, 144, 41, 89, 122, 120, 82, 166, 185, 57, 33, 171, 56, 148, 72, 167, 2, 31, 121, 166, 4, 59, 79, 0, 4, 76, 159, 149, 123, 181, 91, 57, 57, 20, 84, 162, 76, 149, 182, 139, 109, 10, 0, 25, 122, 101, 76, 41, 21, 174, 115, 30, 182, 131, 23, 125, 127, 9, 102, 30, 151, 107, 14, 177, 73, 87, 116, 48, 88, 166, 87, 21, 46, 142, 111, 182, 186, 95, 100, 151, 120, 148, 146, 104, 137, 64, 129, 40, 74, 100, 130, 0, 149, 131, 133, 178, 28, 158, 29, 183, 172, 53, 107, 35, 64, 119, 27, 186, 38, 8, 187, 120, 16, 161, 166, 144, 169, 116, 178, 141, 153, 155, 125, 126, 47, 12, 164, 115, 56, 22, 88, 20, 71, 124, 177, 55, 104, 79, 109, 64, 184, 61, 111, 142, 97, 166, 74, 32, 159, 19, 34, 27, 68, 122, 63, 112, 185, 103, 30, 20, 88, 15, 77, 161, 55, 29, 156, 172, 83, 62, 75, 39, 55, 30, 44, 185, 20, 30, 40, 187, 161, 114, 119, 66, 90, 62, 28, 74, 186, 145, 14, 130, 147, 7, 88, 109, 89, 183, 79, 172, 61, 82, 107, 61, 103, 177, 70, 114, 162, 154, 14, 149, 180, 64, 160, 74, 185, 25, 105, 58, 120, 17, 150, 179, 171, 185, 89, 63, 90, 18, 35, 64, 88, 22, 149, 86, 7, 115, 85, 23, 130, 157, 112, 175, 114, 106, 187, 49, 28, 61, 140, 0, 154};
		vector<int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0])));
		string Arg3 = "Chores";
		verify_case(0, Arg3, canWalkExactly(Arg0, Arg1, Arg2));
	}
	void test_case_1() { int Arg0 = 3; int Arr1[] = {0, 1, 1, 2}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {1, 0, 2, 0}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); string Arg3 = "Freedom"; verify_case(1, Arg3, canWalkExactly(Arg0, Arg1, Arg2)); }
	void test_case_2() { int Arg0 = 4; int Arr1[] = {0, 2, 2, 3, 0}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {2, 0, 3, 0, 1}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); string Arg3 = "Chores"; verify_case(2, Arg3, canWalkExactly(Arg0, Arg1, Arg2)); }
	void test_case_3() { int Arg0 = 10; int Arr1[] = {0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 9, 9}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {4, 5, 6, 7, 9, 4, 6, 0, 1, 3, 8, 4, 6, 1, 4, 8, 1, 7, 8, 1, 4, 2, 5, 6, 8}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); string Arg3 = "Chores"; verify_case(3, Arg3, canWalkExactly(Arg0, Arg1, Arg2)); }
	void test_case_4() { int Arg0 = 8; int Arr1[] = {0, 1, 4, 6, 7, 5, 2, 3, 0}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {1, 4, 6, 7, 5, 2, 3, 0, 7}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); string Arg3 = "Freedom"; verify_case(4, Arg3, canWalkExactly(Arg0, Arg1, Arg2)); }
	void test_case_5() { int Arg0 = 2000; int Arr1[] = {}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); string Arg3 = "Chores"; verify_case(5, Arg3, canWalkExactly(Arg0, Arg1, Arg2)); }

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	WalkingToSchool ___test;
	___test.run_test(0);

	return 0;
}
// END CUT HERE


