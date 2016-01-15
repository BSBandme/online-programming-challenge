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

template <class T> inline T abs1(T a) {return a < 0 ? -a : a;}


inline int jud(double a, double b){
	if(abs(a) < eps && abs(b) < eps) return 0;
	else if(abs1(a - b) / abs1(a) < eps) return 0;
	if(a < b) return -1;
	return 1;
}

//.........................mi.......feng......xian.......xia.......jin.......zhi.......hack...............................................


struct WaterTank{
	double minOutputRate(vector <int> t, vector <int> x, int C){
		double be = 0, en = 1000000;
		int n = t.size();
		while(fabs(be - en) > 1e-6) {
			double rc = 0;
			double mid = (be + en) / 2;
			bool flag = 1;
			for(int i = 0; i < n; i++) {
				rc -= t[i] * (mid - x[i]);
				rc = max(rc, 0.0);
				if(jud(rc, C * 1.0) > 0) {
					flag = 0;
					break;
				}
			}
			if(flag) en = mid;
			else be = mid;
		}
		return (be + en) / 2;
	}
	


};




// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor


// Powered by FileEdit
// Powered by TZTester 1.01 [25-Feb-2003]
// Powered by CodeProcessor

