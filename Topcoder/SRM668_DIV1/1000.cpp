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

#define DEBUG
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

/*
 *   dp[i][j] : number of ways that length is i and convert memory from 0 to j;
 *   rdp[i][j] : number of ways that length is i and convert memory from j to 0;
 *   kdp[i][j] : number of ways that length is i and convert memory from j to 0, with "[]" surrounded;
 */

const int maxn = 150;
int cdp[maxn][maxn * 2], cdp1[maxn][maxn * 2], cdp2[maxn][maxn * 2], cdp3[maxn][maxn * 2];
int *dp[maxn], *rdp[maxn], *kdp[maxn], *chun[maxn];

struct Brainstuck{
	int countPrograms(int n, int m){
		memset(cdp, 0, sizeof(cdp));
		memset(cdp1, 0, sizeof(cdp1));
		memset(cdp2, 0, sizeof(cdp1));
		memset(cdp3, 0, sizeof(chun));
		for(int i = 0; i < maxn; i++) {
			dp[i] = cdp[i] + maxn;
			rdp[i] = cdp1[i] + maxn;
			kdp[i] = cdp2[i] + maxn;
			chun[i] = cdp3[i] + maxn;
		}
		dp[0][0] = 1;
		kdp[0][0] = 1;
		rdp[0][0] = 1;
		for(int i = 1; i <= n; i++) {
			add(kdp[i][0], 1ll * kdp[i - 1][0] * 2 % m, m);
			for(int j = 0; j + 2 <= i; j++)
				add(kdp[i][0], 1ll * kdp[i - j - 2][0] * kdp[j][0] % m, m);
		}

		chun[0][0] = 1;
		for(int i = 1; i <= n; i++) {
			for(int j = -n; j <= n; j++) {
				add(dp[i][j], dp[i - 1][j - 1], m);
				add(dp[i][j], dp[i - 1][j + 1], m);
				add(rdp[i][j], rdp[i - 1][j - 1], m);
				add(rdp[i][j], rdp[i - 1][j + 1], m);
				add(chun[i][j], chun[i - 1][j - 1], m);
				add(chun[i][j], chun[i - 1][j + 1], m);
			}

			for(int k = 0; k + 2 <= i; k++) {
				for(int j = -n; j <= n; j++)
					add(dp[i][0], 1ll * kdp[k][j] * dp[i - k - 2][j] % m, m);
			}
			for(int j = -n; j <= n; j++) {
				for(int k = 2; k <= i; k++)
					add(rdp[i][j], 1ll * kdp[k - 2][j] * dp[i - k][0] % m, m);
			}
			for(int j = -n; j <= n; j++) if(j != 0) {
				add(kdp[i][j], rdp[i][j], m);
				for(int k = 2; abs(1ll * k * j) <= n; k++)
					add(kdp[i][k * j], chun[i][j], m);
			}
		}

		int ans = 0;
		for(int i = -n; i <= n; i++)
			add(ans, dp[n][i], m);

		return ans;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); if ((Case == -1) || (Case == 5)) test_case_5(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const int &Expected, const int &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() { int Arg0 = 2; int Arg1 = 1000000000; int Arg2 = 5; verify_case(0, Arg2, countPrograms(Arg0, Arg1)); }
	void test_case_1() { int Arg0 = 3; int Arg1 = 1000000000; int Arg2 = 12; verify_case(1, Arg2, countPrograms(Arg0, Arg1)); }
	void test_case_2() { int Arg0 = 5; int Arg1 = 1000000000; int Arg2 = 92; verify_case(2, Arg2, countPrograms(Arg0, Arg1)); }
	void test_case_3() { int Arg0 = 16; int Arg1 = 1000000000; int Arg2 = 55450070; verify_case(3, Arg2, countPrograms(Arg0, Arg1)); }
	void test_case_4() { int Arg0 = 13; int Arg1 = 163; int Arg2 = 64; verify_case(4, Arg2, countPrograms(Arg0, Arg1)); }
	void test_case_5() { int Arg0 = 100; int Arg1 = 1000000000; int Arg2 = 875563034; verify_case(5, Arg2, countPrograms(Arg0, Arg1)); }

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	Brainstuck ___test;
	___test.run_test(3);

	return 0;
}
// END CUT HERE

