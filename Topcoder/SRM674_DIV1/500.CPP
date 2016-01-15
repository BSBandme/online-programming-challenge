//#define DEBUG
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

const int mod = MOD;

const int maxn = 200100;
pll arr[maxn];
int loc[maxn];
pll que[maxn];
pll l[maxn], r[maxn];
int nl, nr;
int n;
set <int> have;

struct FindingKids{
	long long getSum(int n, int q, int A, int B, int C){
		ll a = A, b = B, c = C;
		::n = n;
		have.clear();
		for(int i = 0; i < n; i++) {
			a = (a * b + c) % mod;
			ll p = a % (mod - n + i + 1);
			if(have.find(p) != have.end())
				p = mod - n + i;
			have.insert(p);
			arr[i] = mpr(p, i);
		}
		for(int i = 0; i < q; i++) {
			a = (a * b + c) % mod;
			que[i].first = a % n;
			a = (a * b + c) % mod;
			que[i].second = a;
		}
		sort(arr, arr + n);
		nl = 0, nr = 0;
		for(int i = 0; i < n; i++) {
			if(arr[i].first % 2 == 0)
				r[nr++] = arr[i];
			else
				l[nl++] = arr[i];
			loc[arr[i].second] = i;
		}

		ll ans = 0;
		for(int i = 0; i < q; i++) {
			int rloc = loc[que[i].first];
			int pl = 0, pr = nl - 1;
			while(pl != pr) {
				int mid = (pl + pr) / 2;
				int rl = 0, rr = nr;
				while(rl < rr) {
					int rmid = (rl + rr) / 2;
					if(l[mid].first - que[i].second * 2 > r[rmid].first)
						rl = rmid + 1;
					else rr = rmid;
				}
				if(mid + rl >= rloc)
					pr = mid;
				else pl = mid + 1;
			}
			int pans = pl;
			int rl = 0, rr = nr;
			while(rl < rr) {
				int rmid = (rl + rr) / 2;
				if(l[pl].first - que[i].second * 2 > r[rmid].first)
					rl = rmid + 1;
				else rr = rmid;
			}
			pans += rl;
			if(pans == rloc) {
				ans += abs(l[pl].first - que[i].second);
				continue;
			}
			pl = 0, pr = nr - 1;
			while(pl != pr) {
				int mid = (pl + pr) / 2;
				int rl = 0, rr = nl;
				while(rl < rr) {
					int rmid = (rl + rr) / 2;
					if(r[mid].first + que[i].second * 2 > l[rmid].first)
						rl = rmid + 1;
					else rr = rmid;
				}
				pans = mid + rl;
				if(mid + rl >= rloc)
					pr = mid;
				else pl = mid + 1;
			}
			ans += abs(r[pl].first + que[i].second);
			//
			#ifdef DEBUG //......................................................................................................
			pans = pl;
			rl = 0, rr = nl;
			while(rl < rr) {
				int rmid = (rl + rr) / 2;
				if(r[pl].first + que[i].second * 2 > l[rmid].first)
					rl = rmid + 1;
				else rr = rmid;
			}
			pans += rl;
			if(pans != rloc) {
				cerr << 1 << endl;
			}
			#endif //...........................................................................................................

		}

		return ans;
	}
	


};

// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor


// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

