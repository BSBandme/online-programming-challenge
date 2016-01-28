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


typedef vector <int> vi;

const int maxn = 1010;
int fa[maxn], arr[maxn], larr;
int used[maxn][maxn];

struct BearSpans{
	vector <int> findAnyGraph(int n, int m, int k){
		vi ans;
		if(m < n - 1) {
			return vi(1, -1);
		}
		memset(used, 0, sizeof(used));
		memset(fa, -1, sizeof(fa));
		for(int i = 0; i < k - 1; i++) {
			int no = -1;
			for(int j = 0; j < n; j++) if(fa[j] == -1){
				if(no == -1) no = j;
				else  {
					ans.push_back(no + 1);
					ans.push_back(j + 1);
					used[no][j] = used[j][no] = 1;
					fa[j] = no;
					no = -1;
				}
			}
			if(no != -1) for(int j = 0; j < n; j++) if(fa[j] == -1) {
				ans.push_back(j + 1);
				ans.push_back(no + 1);
				used[j][no] = used[no][j] = 1;
				fa[no] = j;
				no = -1;
				break;
			}
		}
		int mu = -1, rcnt = 0;
		for(int i = 0; i < n; i++) if(fa[i] == -1) {
			mu = i;
			rcnt++;
		}
		if(rcnt <= 1)
			return vi(1, -1);
		for(int i = 0; i < n; i++) if(fa[i] == -1 && i != mu) {
			ans.push_back(i + 1);
			ans.push_back(mu + 1);
			used[i][mu] = used[mu][i] = 1;
			fa[i] = mu;
		}
		int rm = n - 1;
		for(int i = 0; i < n && rm < m; i++) for(int j = i + 1; j < n && rm < m; j++) if(!used[i][j]) {
			rm++;
			used[i][j] = 1;
			used[j][i] = 1;
			ans.push_back(i + 1);
			ans.push_back(j + 1);
		}
		return  ans ;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); if ((Case == -1) || (Case == 5)) test_case_5(); if ((Case == -1) || (Case == 6)) test_case_6(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const vector <int> &Expected, const vector <int> &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: " << print_array(Expected) << endl; cerr << "\tReceived: " << print_array(Received) << endl; } }
	void test_case_0() { int Arg0 = 17; int Arg1 = 22; int Arg2 = 1; int Arr3[] = {1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9, 1, 10, 1, 11, 1, 12, 1, 13, 1, 14, 1, 15, 1, 16, 1, 17, 2, 3, 2, 8, 3, 9, 8, 9, 10, 12, 12, 14 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(0, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }
	void test_case_1() { int Arg0 = 9; int Arg1 = 22; int Arg2 = 3; int Arr3[] = {6, 5, 7, 6, 1, 2, 3, 4, 8, 9, 3, 5, 7, 4, 1, 9, 6, 2, 6, 1, 1, 8, 4, 5 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(1, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }
	void test_case_2() { int Arg0 = 9; int Arg1 = 12; int Arg2 = 7; int Arr3[] = {-1 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(2, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }
	void test_case_3() { int Arg0 = 1000; int Arg1 = 999; int Arg2 = 970; int Arr3[] = {-1 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(3, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }
	void test_case_4() { int Arg0 = 2; int Arg1 = 1; int Arg2 = 1; int Arr3[] = {1, 2 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(4, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }
	void test_case_5() { int Arg0 = 3; int Arg1 = 2; int Arg2 = 1; int Arr3[] = {1, 2, 1, 3 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(5, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }
	void test_case_6() { int Arg0 = 3; int Arg1 = 3; int Arg2 = 2; int Arr3[] = {-1 }; vector <int> Arg3(Arr3, Arr3 + (sizeof(Arr3) / sizeof(Arr3[0]))); verify_case(6, Arg3, findAnyGraph(Arg0, Arg1, Arg2)); }

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	BearSpans ___test;
	___test.run_test(-1);

	return 0;
}
// END CUT HERE 


