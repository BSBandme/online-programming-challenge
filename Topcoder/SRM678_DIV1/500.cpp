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


typedef long long ll;
typedef pair <ll, ll> pll;

#define  MOD 1000000007


#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

const int mod = MOD;
const int maxn = 100100;
pll arr[maxn];
int n;

struct TheEmpireStrikesBack{
	int find(int ax, int bx, int cx, int ay, int by, int cy, int n, int m){
		::n = n;
		arr[0] = make_pair(ax, ay);
		for(int i = 1; i < n; i++)  {
			arr[i].first = (arr[i - 1].first * bx + cx) % mod;
			arr[i].second = (arr[i - 1].second * by + cy) % mod;
		}
		sort(arr, arr + n);
		int rx = 0;
		for(int i = 0; i < n; i++) {
			while(rx && arr[rx - 1].second <= arr[i].second)
				rx--;
			arr[rx++] = arr[i];
		}
		n = rx;
		int be = 0, en = MOD;
		while(be != en) {
			int mid = (be + en) / 2;
			int cnt = 0;
			for(int i = 0; i < n; ) {
				int ri = i;
				for(; ri < n - 1 && arr[ri + 1].second + mid >= arr[i].second; ri++);
				cnt++;
				i = ri;
				while(i < n && arr[i].first <= arr[ri].first + mid)
					i++;
			}
			if(cnt <= m) en = mid;
			else be = mid + 1;
		}

		return be ;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const int &Expected, const int &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() { int Arg0 = 2; int Arg1 = 2; int Arg2 = 2; int Arg3 = 2; int Arg4 = 2; int Arg5 = 2; int Arg6 = 2; int Arg7 = 1; int Arg8 = 0; verify_case(0, Arg8, find(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)); }
	void test_case_1() { int Arg0 = 2; int Arg1 = 2; int Arg2 = 2; int Arg3 = 2; int Arg4 = 4; int Arg5 = 1000000000; int Arg6 = 2; int Arg7 = 1; int Arg8 = 1; verify_case(1, Arg8, find(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)); }
	void test_case_2() { int Arg0 = 1; int Arg1 = 3; int Arg2 = 8; int Arg3 = 10000; int Arg4 = 10; int Arg5 = 999910000; int Arg6 = 3; int Arg7 = 1; int Arg8 = 30; verify_case(2, Arg8, find(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)); }
	void test_case_3() { int Arg0 = 0; int Arg1 = 0; int Arg2 = 0; int Arg3 = 0; int Arg4 = 0; int Arg5 = 0; int Arg6 = 100000; int Arg7 = 1000; int Arg8 = 0; verify_case(3, Arg8, find(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)); }
	void test_case_4() { int Arg0 = 10; int Arg1 = 20; int Arg2 = 30; int Arg3 = 40; int Arg4 = 50; int Arg5 = 60; int Arg6 = 100000; int Arg7 = 10; int Arg8 = 15720; verify_case(4, Arg8, find(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)); }

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	TheEmpireStrikesBack ___test;
	___test.run_test(-1);

	return 0;
}
// END CUT HERE 


