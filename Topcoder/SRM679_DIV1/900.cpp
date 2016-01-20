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

const int maxn = 550;
ll cnt[maxn][maxn], sum[maxn][maxn];


struct BagAndCards{
	int getHash(int n, int m, int x, int a, int b, int c, string is){
		memset(cnt, 0, sizeof(cnt));
		memset(sum, 0, sizeof(sum));
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j ++ ) {
				cnt[i][j] = x;
				for(int k = 0; k < m; k++)
					if(is[j + k] == 'Y')
						add(sum[i][k], x);
				x = ((1ll * x * a + b) ^ (1ll * c)) % MOD;
			}
		}

		ll ans = 0;
		for(int i = 0; i < n; i++) for(int j = i + 1; j < n; j++) {
			ll rsum = 0;
			for(int k = 0; k < m; k++)
				add(rsum, cnt[i][k] * sum[j][k] % MOD);
			ans ^= rsum;
		}
		return ans;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const int &Expected, const int &Received) {
		cerr << "Test Case #" << Case << "...";
		if (Expected == Received)
			cerr << "PASSED" << endl;
		else {
			cerr << "FAILED" << endl;
			cerr << "\tExpected: \"" << Expected << '\"' << endl;
			cerr << "\tReceived: \"" << Received << '\"' << endl;
		}
	}
	void test_case_0() {
		int Arg0 = 2;
		int Arg1 = 4;
		int Arg2 = 1;
		int Arg3 = 1;
		int Arg4 = 0;
		int Arg5 = 0;
		string Arg6 = "NNYYNYN";
		int Arg7 = 9;
		verify_case(0, Arg7, getHash(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6));
	}
	void test_case_1() {
		int Arg0 = 3;
		int Arg1 = 5;
		int Arg2 = 1;
		int Arg3 = 1;
		int Arg4 = 1;
		int Arg5 = 2;
		string Arg6 = "NNYYNYNYN";
		int Arg7 = 1532;
		verify_case(1, Arg7, getHash(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6));
	}
	void test_case_2() {
		int Arg0 = 10;
		int Arg1 = 20;
		int Arg2 = 111;
		int Arg3 = 222;
		int Arg4 = 333;
		int Arg5 = 444;
		string Arg6 = "NNNNNYYYNNNYYYYYYNNYYYYNNNNYNNYYYNNNYYN";
		int Arg7 = 450750683;
		verify_case(2, Arg7, getHash(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6));
	}
	void test_case_3() {
		int Arg0 = 2;
		int Arg1 = 2;
		int Arg2 = 1;
		int Arg3 = 1;
		int Arg4 = 0;
		int Arg5 = 0;
		string Arg6 = "NNY";
		int Arg7 = 1;
		verify_case(3, Arg7, getHash(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6));
	}

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	BagAndCards ___test;
	___test.run_test(-1);

	return 0;
}
// END CUT HERE 


