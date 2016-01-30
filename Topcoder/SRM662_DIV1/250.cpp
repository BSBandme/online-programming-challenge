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
pii arr[maxn];

struct FoxesOfTheRoundTable{
  vector <int> minimalDifference(vector <int> h){
    int n = h.size();
    for(int i = 0; i < n; i++)
      arr[i] = mpr(h[i], i);
    sort(arr, arr + n);

    vi ans, ans1;
    for(int i = 0; i < n; i += 2) {
      ans.push_back(arr[i].second);
    }
    for(int i = 1; i < n; i += 2)
      ans1.push_back(arr[i].second);
    for(int i = ans1.size() - 1; i >= 0; i--)
      ans.push_back(ans1[i]);
    return vector <int>(ans) ;
  }
  


};

