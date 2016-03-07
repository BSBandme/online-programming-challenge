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


const int maxn = 55;
struct edge {
	int to, nxt;
} e[maxn * 2];
int head[maxn], le;
int cnt[maxn];
int st[maxn], lst, in[maxn], iscir[maxn];
int du[maxn];

vi cir;
int n;

void addedge(int a, int b) {
	e[le].to = b;
	e[le].nxt = head[a];
	head[a] = le++;
}
void dfs(int no, int fae) {
	if(in[no]) {
		cir.clear();
		for(int j = lst - 1; j >= 0; j--) {
			iscir[st[j]] = 1;
			if(st[j] == no) break;
		}
		return;
	}
	in[no] = 1;
	st[lst++] = no;
	for(int i = head[no]; i != -1; i = e[i].nxt) if(i != fae) {
		dfs(e[i].to, i ^ 1);
	}
	lst--;
	in[no] = 0;
}
int dfs1(int no, int fae) {
	int ans = 0;
	for(int i = head[no]; i != -1; i = e[i].nxt) if(i != fae && !iscir[e[i].to]) {
		ans += dfs1(e[i].to, i ^ 1);
	}
	if(ans == 0 && !iscir[no]) ans = 1;
	return ans;
}

struct SuccessfulMerger{
	int minimumMergers(vector <int> road){
		n = road.size();
		memset(head, -1, sizeof(head));
		le = 0;
		memset(du, 0, sizeof(du));
		memset(iscir, 0, sizeof(iscir));
		for(int i = 0; i < n; i++) {
			du[i]++;
			du[road[i]]++;
			addedge(i, road[i]);
			addedge(road[i], i);
		}
		int ans = 0;
		for(int i = 0; i < n; i++) if(du[i] == 1)
			ans++;
		ans = n - (ans + 1);

		dfs(0, -1);
		for(int i = 0; i < n; i++) if(iscir[i]){
			if(dfs1(i, -1) == 0) {
				ans--;
				break;
			}
		}

		return ans;
	}
	
// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); if ((Case == -1) || (Case == 5)) test_case_5(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const int &Expected, const int &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() { int Arr0[] = {2, 2, 1, 1, 1}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arg1 = 1; verify_case(0, Arg1, minimumMergers(Arg0)); }
	void test_case_1() { int Arr0[] = {3, 4, 5, 4, 5, 3}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arg1 = 2; verify_case(1, Arg1, minimumMergers(Arg0)); }
	void test_case_2() { int Arr0[] = {1, 0, 1, 0, 0, 0, 1, 0, 1, 1}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arg1 = 1; verify_case(2, Arg1, minimumMergers(Arg0)); }
	void test_case_3() { int Arr0[] = {1, 2, 3, 0}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arg1 = 2; verify_case(3, Arg1, minimumMergers(Arg0)); }
	void test_case_4() { int Arr0[] = {31, 19, 0, 15, 30, 32, 15, 24, 0, 20, 40, 1, 22, 21, 39, 28, 0, 23, 15, 5, 19, 22, 6, 32, 8, 38, 35, 30, 4, 28, 32, 18, 18, 9, 22, 41, 20, 18, 6, 25, 41, 34, 4}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); int Arg1 = 25; verify_case(4, Arg1, minimumMergers(Arg0)); }
	void test_case_5() {
		int Arr0[] = {9, 3, 7, 8, 12, 18, 2, 11, 19, 7, 13, 13, 1, 19, 19, 4, 10, 19, 6, 3};
		vector<int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
		int Arg1 = 12;
		verify_case(5, Arg1, minimumMergers(Arg0));
	}

// END CUT HERE


};

// BEGIN CUT HERE
int main(){
	SuccessfulMerger ___test;
	___test.run_test(5);

	return 0;
}
// END CUT HERE


