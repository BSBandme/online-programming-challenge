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

const int maxn = 3010;

struct ANewHope{
	int count(vector <int> st, vector <int> en, int d){
		int n = st.size();
		int rd = n - d;
		int delta[maxn];
		memset(delta, 0, sizeof(delta));
		for(int i = 0; i < n; i++) {
			delta[st[i]] += i;
			delta[en[i]] -= i;
		}
		int maxv = 0;
		for(int j = 1; j <= n; j++)
			maxv = max(delta[j], maxv);
		if(maxv == 0) return 1;

		return (maxv - 1) / rd + 2;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

