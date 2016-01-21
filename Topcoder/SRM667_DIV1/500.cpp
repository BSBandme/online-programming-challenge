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


const double pi = acos(0.0) * 2.0;
const double eps = 1e-12;
const int step[8][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

template <class T> inline T abs1(T a) {return a < 0 ? -a : a;}

inline int jud(double a, double b){
	if(abs1(a - b) / abs1(a) < eps) return 0;
	if(a < b) return -1;
	return 1;
}



#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

long double calc(double p, double l) {
	if(l == 0) return p;
	long double xi = p;
	for(int i = 0; i < l - 1; i++) {
		xi = 1 - xi * (1 - p);
		xi = p / xi;
	}
	return xi;
}

const int maxn = 2010;
double dp[maxn][maxn][2];

struct CatsOnTheCircle{
	double getProb(int n, int k, int p){
		if(k == 0) {
			if(n == 1) return 1;
			return 0;
		}
		if(n == 2) {
			if(k == 2) return 1;
			else return 0;
		}
		double rp = 1.0 * p / 1000000000;
		memset(dp, 0, sizeof(dp));
		dp[2][0][0] = 1 - rp;
		dp[2][1][1] = rp;
		for(int i = 2; i < min(n - 1, 2000); i++) {
			double p1 = calc(rp, i);
			double p2 = calc(1 - rp, i);
			for(int j = 0; j < i; j++) {
				dp[i + 1][j][0] += dp[i][j][0] * p2;
				dp[i + 1][j + 1][1] += dp[i][j][0] * (1 - p2);
				dp[i + 1][j + 1][1] += dp[i][j][1] * p1;
				dp[i + 1][j][0] += dp[i][j][1] * (1 - p1);
			}
		}
		if(n - 1 <= 2000) {
			return dp[n - 1][k - 1][0] + dp[n - 1][k - 1][1];
		}
		if(k > 2000 && n - k > 2000)
			return 0;
		if(n - k <= 2000) {
			if((dp[2000][2000 - (n - k)][0] + dp[2000][2000 - (n - k)][1]) == 0)
				return (dp[2000][2000 - (n - k)][0] + dp[2000][2000 - (n - k)][1]);
			return (dp[2000][2000 - (n - k)][0] + dp[2000][2000 - (n - k)][1]) / pow((dp[1000][1000 - (n - k)][0] + dp[1000][1000 - (n - k)][1]) / (dp[2000][2000 - (n - k)][0] + dp[2000][2000 - (n - k)][1]), (n - 2000) / 1000);
		}
		else {
			if(dp[2000][k - 1][0] + dp[2000][k - 1][1] == 0)
				return (dp[2000][k - 1][0] + dp[2000][k - 1][1]);
			return (dp[2000][k - 1][0] + dp[2000][k - 1][1]) / pow((dp[1000][k - 1][0] + dp[1000][k - 1][1]) / (dp[2000][k - 1][0] + dp[2000][k - 1][1]), (n - 2000) / 1000);
		}
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

