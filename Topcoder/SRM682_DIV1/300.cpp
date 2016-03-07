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


#define debug
//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................

const int maxn = 55;
pair <pii, int> arr[maxn];
pair <pii, int> rarr[maxn];
ll rcnt[maxn];
int n, m;

struct FleetFunding{
	int maxShips(int m, vector <int> k, vector <int> a, vector <int> b){
		::m = m;
		n = k.size();
		for(int i = 0; i < n; i++) {
			arr[i] = mpr( mpr(a[i] - 1, b[i] - 1), k[i] );
		}
		sort(arr, arr + n);
		int be = 0, en = 1000000000;
		while(be != en) {
			for(int i = 0; i < n; i++)
				rcnt[i] = arr[i].second;
			int mid = (be + en + 1) / 2;
			set <pii> s;
			bool flag = 1;
			for(int i = 0, link = 0; i < m; i++) {
				while(s.size() && s.begin()->first < i)
					s.erase(s.begin());
				for(; link < n && arr[link].first.first <= i; link++)
					s.insert(mpr(arr[link].first.second, link));
				int left = mid;
				for(; s.size() && left; ) {
					auto no = s.begin()->second;
					s.erase(s.begin());
					if(rcnt[no] >= left) {
						rcnt[no] -= left;
						left = 0;
						if(rcnt[no])
							s.insert(mpr(arr[no].first.second, no));
						break;
					} else
						left -= rcnt[no];
				}
				if(left) {
					flag = 0;
					break;
				}
			}
			if(!flag) en = mid - 1;
			else be = mid;
		}

		return be;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor


// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

