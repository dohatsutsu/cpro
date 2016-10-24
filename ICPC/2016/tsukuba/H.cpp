#include<bits/stdc++.h>
using namespace std;
#define MAX_N 100005
#define MAX_M 100005
 
int n,m,ans;
int x[MAX_M],y[MAX_M],w[MAX_M];
vector<int> G[MAX_N];
vector<int> g[MAX_N];
vector<int> group[MAX_N];
int id[MAX_N];
bool visited[MAX_N];
vector<int> vd;
int dp[MAX_N];
int d[MAX_N];
 
bool dfs(int pos,int prev,int root){
  id[pos]=root;
  visited[pos]=true;
  group[root].push_back(pos);
  for(int i=0;i<(int)G[pos].size();i++){
    int to=G[pos][i];
    if(to==prev)continue;
    if(i&&to==G[pos][i-1])return true;
    if(visited[to])return true;
    if(dfs(to,pos,root))return true;
  }
  return false;
}
 
bool check(){
  memset(visited,false,sizeof(visited));
  for(int i=0;i<n;i++)
    if(!visited[i]&&dfs(i,-1,i))return true;
 
  vector<int> deg(n);
  for(int pos=0;pos<n;pos++)
    for(int to:g[pos])
      deg[ id[to] ]++;
 
  queue<int> Q;
  for(int i=0;i<n;i++)
    if(deg[i]==0&&id[i]==i)
      Q.push(i);
   
  while(!Q.empty()){
    int root=Q.front();Q.pop();
    vd.push_back(root);
    for(int pos:group[root]){
      for(int to:g[pos]){
        deg[ id[to] ]--;
        if(deg[ id[to] ]==0)Q.push(id[to]);
      }
    }
  }
  for(int i=0;i<n;i++)
    if(deg[i]>0)return true;
  return false;
}
 
int solve(int S){
  for(int pos:group[id[S]])d[pos]=-1;
  queue<int> Q;
  Q.push(S);
  d[S]=0;
  while(!Q.empty()){
    int pos=Q.front();Q.pop();
    for(int to:G[pos]){
      if(d[to]!=-1)continue;
      d[to]=d[pos]+1;
      Q.push(to);
    }
  }
  int res=S;
  for(int pos:group[id[S]]){
    if(d[pos]+dp[pos]>d[res]+dp[res])res=pos;
    ans=max(ans,d[pos]+dp[S]);
    for(int to:g[pos])dp[to]=max(dp[to],d[pos]+dp[S]+1);
  }
  return res;
}
 
int main(){
  scanf("%d %d",&n,&m);
  for(int i=0;i<m;i++){
    scanf("%d %d %d",&x[i],&y[i],&w[i]);
    x[i]--,y[i]--;
    if(w[i]==1){
      g[x[i]].push_back(y[i]);
    }
    if(w[i]==2){
      G[x[i]].push_back(y[i]);
      G[y[i]].push_back(x[i]);
    }
  }
  if(check()){
    printf("Infinite\n");
    return 0;
  }
  for(int root : vd )
    if(root==id[root])
      solve(solve(solve(root)));
  printf("%d\n",ans);  
  return 0;
}
