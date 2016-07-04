#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef complex<ll> P;
  
namespace std{
  bool operator < (const P &a,const P &b){
    return (a.real()==b.real()?a.imag()<b.imag():a.real()<b.real());
  }
};
  
int ccw(P a,P b,P c){
  b-=a;c-=a;
  P k=c*conj(b);
  if(k.imag() > 0)return 1;
  if(k.imag() < 0)return -1;
  if(k.real() < 0)return 2;
  if(norm(b) > norm(c))return -2;
  return 0;
}
  
double dist(P a,P b){
  double x=real(a-b);
  double y=imag(a-b);
  return sqrt(x*x+y*y);
}
  
int n,m,size;
P s[200005],t[200005],u[200005];
vector<int> G[200005];
map< P , bool > mp;
 
double solve(){
  double res=0,sum=0;
  int k=0;
  for(int i=0;i<size;i++){
    if(i&&t[i]==t[i-1])continue;
    while(k>=2 && ccw(u[k-2],u[k-1],t[i])!=1){
      sum-=dist(u[k-2],u[k-1]);
      k--;
    }
    if(k>=1)sum+=dist(u[k-1],t[i]);      
    u[k++]=t[i];
    if(mp[ t[i] ]) res+=sum;
  }
  return res;
}
 
int main(){
  int Tc,x,y;
  scanf("%d",&Tc);
  for(int tc=1;tc<=Tc;tc++){
    scanf("%d",&n);
    mp.clear();
    for(int i=0;i<n;i++)G[i].clear();
     
    for(int i=0;i<n;i++){
      scanf("%d %d",&x,&y);
      s[i]=P(x,y);
    }
     
    scanf("%d",&m);
    for(int i=0;i<m;i++){
      scanf("%d %d",&x,&y);
      G[x].push_back(y);
    }
     
    size=0;
    for(int i=0;i<n;i++){
      t[size++]=s[i];
 
      P A=s[i],B=s[(i+1)%n],AB=(B-A)/abs(B-A);
      sort(G[i].begin(),G[i].end());
      for(int j=0;j<(int)G[i].size();j++){
        ll k=G[i][j];
        mp[A+AB*k]=true;
        t[size++]=A+AB*k;
      }
    }
    reverse(t+1,t+size);
    printf("%.1f\n",solve());
  
  }
  return 0;
}
