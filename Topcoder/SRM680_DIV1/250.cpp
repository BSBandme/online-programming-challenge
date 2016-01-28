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


const int maxn = 1010;
int dp[maxn][55][55];
int cons[maxn];

struct BearFair{
	string isFair(int n, int b, vector <int> u, vector <int> q){
		if(n & 1)
			return "unfair";
		memset(cons, -1, sizeof(cons));
		for(int i = 0; i < (int)u.size(); i++) {
			if(cons[u[i]] != -1 && cons[u[i]] != q[i])
				return "unfair";
			cons[u[i]] = q[i];
		}
		if(cons[b] != -1)
			if(cons[b] != n)
				return "unfair";
		memset(dp, 0, sizeof(dp));
		dp[0][0][0] = 1;
		for(int i = 0; i < b; i++) {
			int bej = 0, enj = min(i, n);
			if(cons[i] != -1)
				bej = cons[i], enj = cons[i];
			for(int j = bej; j <= enj; j++) for(int k = 0; k <= j; k++) if(dp[i][j][k]){
				dp[i + 1][j + 1][k + (i & 1)] = 1;
				dp[i + 1][j][k] = 1;
			}
		}
		if(dp[b][n][n / 2])
			return "fair";
		return "unfair";
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

