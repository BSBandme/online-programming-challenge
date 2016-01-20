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

const double pi = acos(0.0) * 2.0;
const double eps = 1e-10;

template <class T> inline T abs1(T a) {return a < 0 ? -a : a;}

inline int jud(double a, double b){
	if(abs(a) < eps && abs(b) < eps) return 0;
	else if(abs1(a - b) / abs1(a) < eps) return 0;
	if(a < b) return -1;
	return 1;
}


#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................


struct point {
	union {
		double co[2];
		struct {
			double x, y;
		};
	};
	double ang;
	point(double a = 0, double b = 0) : x(a), y(b) {ang = 0;}

	point operator + (const point a) {
		point ans;
		ans.x = x + a.x;
		ans.y = y + a.y;
		return ans;
	}
	point operator - (const point a) {
		point ans;
		ans.x = x - a.x;
		ans.y = y - a.y;
		return ans;
	}
	double operator * (const point a) const {
		return x * a.x + y * a.y;
	}
	double operator % (const point a) {
		return x * a.y - y * a.x;
	}
	point operator * (double p) {
		point ans;
		ans.x = x * p;
		ans.y = y * p;
		return ans;
	}

};
double getang(double x, double y) {
	if(x == 0) {
		if(y > 0) return pi / 2;
		else return pi / 2 * 3;
	}
	if(x > 0) {
		if(y < 0) return atan(y / x) + pi * 2;
		else return atan(y / x);
	} else {
		return atan(y / x) + pi;
	}
}
double getang(point p) {
	return getang(p.x, p.y);
}
bool in(point cen, point a, point b) {
	return jud((a - cen) % (b - cen), 0) >= 0;
}

//判断两线段是否相交
bool inter(point a, point b, point c, point d){
	if ( min(a.x, b.x) > max(c.x, d.x) ||
	min(a.y, b.y) > max(c.y, d.y) ||
	min(c.x, d.x) > max(a.x, b.x) ||
	min(c.y, d.y) > max(a.y, b.y) ) return 0;
	double h, i, j, k;
	h = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	i = (b.x - a.x) * (d.y - a.y) - (b.y - a.y) * (d.x - a.x);
	j = (d.x - c.x) * (a.y - c.y) - (d.y - c.y) * (a.x - c.x);
	k = (d.x - c.x) * (b.y - c.y) - (d.y - c.y) * (b.x - c.x);
	return h * i <= eps && j * k <= eps;
}
struct line {
	point a, b;
	double ang;
	line(point a1 = point(0, 0), point b1 = point(1, 0)): a(a1), b(b1) {
		point temp = b1 - a1;
		ang = getang(temp.x, temp.y);
	}
	friend bool contain(line l, point a);
	bool operator == (const line &rline) const {
		return contain(*this, rline.a) && contain(*this, rline.b);
	}
};
inline bool contain(line l, point a) {
	return jud((a - l.a) % (l.b - l.a), 0) == 0;
}
point getinter(line la, line lb) {
	double sa = (lb.a - la.a) % (lb.a - la.b);
	double sb = (lb.b - la.b) % (lb.b - la.a);
	point ans = (lb.a * sb + lb.b * sa) *(1.0 / (sb + sa));
	return ans;
}

const int maxn = 55;
point arrb[maxn], arrr[maxn];
int lb, lr;

bool cmp(const point &a, const point &b) {
	if(a.x == b.x) return a.y < b.y;
	return a.x < b.x;
}

int dp[maxn][maxn];
int can[maxn][maxn];
int sarr[maxn], lsarr;
int base;

bool cmp1(const int &a, const int &b) {
	if(a == base) return a < b;
	if((arrb[a] - arrb[base]) % (arrb[b] - arrb[base]) == 0) {
		return arrb[a].x - arrb[base].x <= arrb[b].x - arrb[base].x &&
				arrb[a].y - arrb[base].y <= arrb[b].y - arrb[base].y;
	}
	return (arrb[a] - arrb[base]) % (arrb[b] - arrb[base]) < 0;
}

struct RedAndBluePoints{
	int find(vector <int> bluex, vector <int> bluey, vector <int> redx, vector <int> redy){
		lb = bluex.size();
		lr = redx.size();
		for(int i = 0; i < lb; i++)
			arrb[i] = point(bluex[i], bluey[i]);
		for(int i = 0; i < lr; i++)
			arrr[i] = point(redx[i], redy[i]);

		sort(arrb, arrb + lb, cmp);

		int ans = 1;
		for(int i = 0; i < lb; i++) {
			lsarr = 0;
			base = i;
			for(int j = i; j < lb; j++) {
				sarr[lsarr++] = j;
			}
			sort(sarr, sarr + lsarr, cmp1);
			for(int j = 0; j < lsarr; j++) for(int k = j + 1; k < lsarr; k++){
				int noj = sarr[j], nok = sarr[k];
				bool flag = 1;
				for(int t = 0; t < lr && flag; t++) {
					if(in(arrr[t], arrb[base], arrb[noj]) && in(arrr[t], arrb[noj], arrb[nok]) && in(arrr[t], arrb[nok], arrb[base]))
						flag = 0;
					if(in(arrr[t], arrb[base], arrb[nok]) && in(arrr[t], arrb[nok], arrb[noj]) && in(arrr[t], arrb[noj], arrb[base]))
						flag = 0;
				}
				if(flag) {
					can[j][k] = 1;
					for(int t = j + 1; t < k; t++) {
						if(in(arrb[sarr[t]], arrb[base], arrb[noj]) && in(arrb[sarr[t]], arrb[noj], arrb[nok]) && in(arrb[sarr[t]], arrb[nok], arrb[base]))
							can[j][k]++;
						else if(in(arrb[sarr[t]], arrb[base], arrb[nok]) && in(arrb[sarr[t]], arrb[nok], arrb[noj]) && in(arrb[sarr[t]], arrb[noj], arrb[base]))
							can[j][k]++;
					}
				} else
					can[j][k] = -1;
			}
			memset(dp, -1, sizeof(dp));
			for(int j = 1; j < lsarr; j++) if(can[0][j] != -1)
				dp[0][j] = 2;
			for(int j = 0; j < lsarr; j++) for(int k = j + 1; k < lsarr; k++) if(dp[j][k] != -1){
				for(int t = k + 1; t < lsarr; t++) if(can[k][t] != -1 && in(arrb[sarr[k]], arrb[sarr[j]], arrb[sarr[t]]))
					dp[k][t] = max(dp[j][k] + can[k][t], dp[k][t]);
			}
			for(int i = 0; i < lsarr; i++) for(int j = i + 1; j < lsarr; j++)
				ans = max(ans, dp[i][j]);
		}

		return ans;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const int &Expected, const int &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() { int Arr0[] = {0,0,10,10}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arr1[] = {0,10,0,10}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {100}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); int Arr3[] = {120}; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); int Arg4 = 4; verify_case(0, Arg4, find(Arg0, Arg1, Arg2, Arg3)); }
	void test_case_1() { int Arr0[] = {0,0,10,10}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arr1[] = {0,10,0,10}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {3}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); int Arr3[] = {4}; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); int Arg4 = 3; verify_case(1, Arg4, find(Arg0, Arg1, Arg2, Arg3)); }
	void test_case_2() { int Arr0[] = {0,0,10,10}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arr1[] = {0,10,0,10}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {3,6}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); int Arr3[] = {2,7}; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); int Arg4 = 2; verify_case(2, Arg4, find(Arg0, Arg1, Arg2, Arg3)); }
	void test_case_3() {
		int Arr0[] = { 0, 0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6 };
		vector<int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
		int Arr1[] = { 0, 2, 4, 6, 0, 2, 4, 6, 0, 2, 4, 6, 0, 2, 4, 6 };
		vector<int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0])));
		int Arr2[] = { 1, 1, 5, 5 };
		vector<int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0])));
		int Arr3[] = { 1, 5, 1, 5 };
		vector<int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0])));
		int Arg4 = 12;
		verify_case(3, Arg4, find(Arg0, Arg1, Arg2, Arg3));
	}
	void test_case_4() { int Arr0[] = {5, 6, 6}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arr1[] = {9, 0, 5}; vector <int> Arg1(Arr1, Arr1 + (sizeof(Arr1) / sizeof(Arr1[0]))); int Arr2[] = {7}; vector <int> Arg2(Arr2, Arr2 + (sizeof(Arr2) / sizeof(Arr2[0]))); int Arr3[] = {6}; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); int Arg4 = 3; verify_case(4, Arg4, find(Arg0, Arg1, Arg2, Arg3)); }

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	RedAndBluePoints ___test;
	___test.run_test(-1);

	return 0;
}
// END CUT HERE 


