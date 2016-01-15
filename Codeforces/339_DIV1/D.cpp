#include <bits/stdc++.h>

using namespace std;

#define mpr make_pair
typedef pair <int, int> pii;

const int maxn = 100100;
int fa[maxn][22], deep[maxn], n;

inline void initfa() {
	for(int i = 1; i <= 20; i++) {
		for(int j = 0; j < n; j++) {
			if(fa[j][i - 1] == -1) fa[j][i] = -1;
			else fa[j][i] = fa[fa[j][i - 1]][i - 1];
		}
	}
}
inline int getfa(int no, int l) {
	for(int i = 0; l; l >>= 1, i++)
		if(l & 1) no = fa[no][i];
	return no;
}
inline int getlca(int a, int b) {
	if(deep[a] > deep[b]) swap(a, b);
	b = getfa(b, deep[b] - deep[a]);
	if(b == a) return a;
	for(int i = 20; i >= 0; i--)
		if(fa[a][i] != fa[b][i]){
			a = fa[a][i];
			b = fa[b][i];
		}
	return fa[a][0];
}
inline int getl(int a, int b) {
	int faa = getlca(a, b);
	return deep[a] + deep[b] - deep[faa] * 2;
}

struct edge {
	int to, nxt;
} e[maxn * 2];
int le, head[maxn];
int sz[maxn];
int arr[maxn], larr;
pii range[maxn];
int nq;
set <pii> s;
set <pair <pii, pii> > pq;
set <int> contain;

void addedge(int a, int b) {
	e[le].to = b;
	e[le].nxt = head[a];
	head[a] = le++;
}

void dfs(int no) {
	sz[no] = 1;
	arr[larr] = no;
	range[no].first = larr++;
	for(int i = head[no]; i != -1; i = e[i].nxt) if(e[i].to != fa[no][0]) {
		fa[e[i].to][0] = no;
		deep[e[i].to] = deep[no] + 1;
		dfs(e[i].to);
	}
	range[no].second = larr - 1;
}

int getnode(int no1, int no2) {
	int lca = getlca(no1, no2);
	if(lca == no1 || lca == no2) {
		if(abs(deep[no1] - deep[no2]) == 1) {
			return -1;
		}
		if(deep[no1] > deep[no2])
			lca = getfa(no1, deep[no1] - deep[no2] - 1);
		else
			lca = getfa(no2, deep[no2] - deep[no1] - 1);
	}
	return lca;
}

int erase(int no) {
	auto it = s.find(mpr(range[no].first, no));
	int l = -1, r = -1;
	if(it != s.begin()) {
		it--;
		l = it->second;
		it++;
	}
	it++;
	if(it != s.end()) {
		r = it->second;
	}
	if(l != -1) {
		int tno = getnode(l, no);
		pq.erase(mpr(mpr(-deep[tno], tno), mpr(l, no)));
	}
	if(r != -1) {
		int tno = getnode(no, r);
		pq.erase(mpr(mpr(-deep[tno], tno), mpr(no, r)));
	}
	if(l != -1 && r != -1) {
		int tno = getnode(l, r);
		if(tno == -1) return 0;
		pq.insert(mpr(mpr(-deep[tno], tno), mpr(l, r)));
	}
	s.erase(mpr(range[no].first, no));
	return 1;
}

int main() {
	memset(head, -1, sizeof(head));
	scanf("%d", &n);
	for(int i = 0; i < n - 1; i++) {
		int a, b;
		scanf("%d%d", &a, &b);
		a--; b--;
		addedge(a, b);
		addedge(b, a);
	}
	fa[0][0] = -1;
	dfs(0);
	initfa();
	scanf("%d", &nq);
	for(int i1 = 0; i1 < nq; i1++) {
		int cnt = 0;
		scanf("%d", &cnt);
		s.clear();
		for(int i = 0; i < cnt; i++) {
			int no;
			scanf("%d", &no);
			no--;
			s.insert(mpr(range[no].first, no));
		}
		contain.clear();
		pq.clear();
		bool flag = 1;
		for(auto it = s.begin(); ;it++) {
			auto it1 = it;
			it1++;
			if(it1 == s.end()) break;
			int no1 = it->second, no2 = it1->second, lca = getnode(no1, no2);
			if(lca == -1) {
				flag = 0;
				break;
			}
			pq.insert(mpr(mpr(-deep[lca], lca), mpr(no1, no2)));
		}
		int ans = 0;
		for(int pcnt = 0; s.size() > 1 && flag; pcnt++) {
			int lca = pq.begin()->first.second;
			set <int> re;
			auto it = pq.begin();
			while(it != pq.end() && it->first.second == lca) {
				re.insert(it->second.first);
				re.insert(it->second.second);
				it++;
			}
			if(contain.find(lca) == contain.end()) {
				contain.insert(lca);
				ans++;
			}
			for(auto it : re) {
				if(deep[it] > deep[lca])
					flag *= erase(it);
				if(!flag) break;
			}
		}
		if(!flag) puts("-1");
		else printf("%d\n", ans);
	}

    return 0;
}


