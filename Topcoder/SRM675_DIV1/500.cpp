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
using namespace std;

#define mpr make_pair
typedef unsigned int ui;
typedef unsigned long long ull;
typedef long long ll;


#define  MOD 1000000007

#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

const int maxn = 100100;
int n, x0, a, b;
int stat[maxn];

struct LimitedMemorySeries1{
	long long getSum(int N, int X0, int A, int B, vector <int> query){
		int bigv = 100000;
		n = N, x0 = X0, a = A, b = B;
		ll ans = 0;
		if(a == 0) {
			if(x0 < b) {
				for(int i = 0; i < (int)query.size(); i++) {
					if(query[i] == 0)
						ans += x0;
					else
						ans += b;
				}
			} else {
				for(int i = 0; i < (int)query.size(); i++) {
					if(query[i] == n - 1)
						ans += x0;
					else
						ans += b;
				}
			}
			return ans;
		}
		if(a == 1 && b == 0) {
			return 1ll * query.size() * x0;
		}

		ll maxv = 0;

		for(ll i = 0, rx = x0; i < n; i++) {
			maxv = max(maxv, rx);
			rx = rx * a + b;
			rx %= MOD;
		}
		bigv = min(1ll * bigv, maxv);

		memset(stat, 0, sizeof(stat));
		int fen = (maxv - 1) / bigv + 1;
		for(ll i = 0, rx = x0; i < n; i++) {
			stat[rx / fen]++;
			rx = rx * a + b;
			rx %= MOD;
		}
		for(int i = 1; i <= bigv; i++)
			stat[i] += stat[i - 1];
		sort(query.begin(), query.end());
		for(int i = 0, j = 0; j < (int)query.size(); i++) if(query[j] < stat[i]) {
			int rj = j + 1;
			int last = 0;
			int ti = i;
			if(i) last = stat[i - 1];
			while(stat[i] - last <= 100000 && rj < (int)query.size()) {
				while(rj < (int)query.size() && query[rj] < stat[i])
					rj++;
				i++;
			}
			int ri = i;
			while(ri && rj && query[rj - 1] < stat[ri - 1])
				ri--;
			vector <int> arr;
			arr.reserve(stat[ri] - last + 10);
			for(ll j = 0, rx = x0; j < n; j++) {
				if(rx / fen >= ti && rx / fen <= ri)
					arr.push_back(rx);
				rx = rx * a + b;
				rx %= MOD;
			}
			sort(arr.begin(), arr.end());
			for(int k = j; k < rj; k++)
				ans += arr[query[k] - last];
			i--;
			j = rj;
		}

		return ans;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

