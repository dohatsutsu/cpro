/*
  問題概要は複雑なので省略
  
  解法　、
  それぞれのアピールを何回行うか2重ループで決めうちする。
　あとはボーカル、ダンス、ビジュアルごとに期待値dpを行う。

*/


#include<bits/stdc++.h>
using namespace std;
#define MAX_N 2005
typedef long double ldouble;

double eps=1e-10;
bool eq(double a,double b){
  return (-eps<a-b&&a-b<eps);
}

int divi(int a,int b){ return (a+b-1)/b; }
 
int n,m;
int a[MAX_N],b[MAX_N],c[MAX_N];
map<int,ldouble> ma,mb,mc;
 
ldouble z[MAX_N];
 
 
ldouble dp[MAX_N][4];
 
ldouble func(int a[MAX_N],map<int,ldouble> &ma,int cnt,ldouble val){
  if(ma.count(cnt))return ma[cnt];
   
  fill( dp[0], dp[MAX_N] , 0.0 );
  dp[1][0]=1;
 
  ldouble w=1.0;
  for(int i=1;i<n;i++){
    int x=divi(a[0]*cnt,a[i]);
     
    ldouble p=(x<=m?z[x]:0.0);
    ldouble q=1.0-p;
    w*=p;
    for(int j=0;j<3;j++){
      if( eq(dp[i][j],0) )continue;
      
      if(!eq(p,0))dp[i+1][j+1]+=dp[i][j]*p;
      if(!eq(q,0))dp[i+1][j]+=dp[i][j]*q;
    }
  }
  ldouble res=(dp[n][0]+dp[n][1]+dp[n][2])*val-w;
  return ma[cnt]=res;  
}
 
ldouble nCr(int n,int r){
  ldouble res=1.0;
  for(int i=0;i<r;i++){
    res/=(ldouble)(i+1);
    res*=(ldouble)(n-i);
  }
  return res;
}
 
int main(){
  scanf("%d %d",&n,&m);
  for(int i=0;i<n;i++){
    scanf("%d %d %d",&a[i],&b[i],&c[i]);
  }
 
  for(int i=0;i<=m;i++){
    z[i]=nCr(m,i);
    for(int j=0;j<i;j++)
      z[i]/=3.0;

    for(int j=0;j<m-i;j++)
      z[i]/=3.0;

    for(int j=0;j<m-i;j++)
      z[i]*=2.0;

  }
  for(int i=m-1;i>=0;i--)z[i]+=z[i+1];
 
  ldouble ans=-1e30;
  for(int i=0;i<=m;i++){
    for(int j=0;j<=m-i;j++){
      int k=m-i-j;
 
      ldouble va=func(a,ma,i,5);
      ldouble vb=func(b,mb,j,3);
      ldouble vc=func(c,mc,k,2);
      ldouble v=va+vb+vc;
 
      ans=max(ans, v);
    }
  }
 
  printf("%.16f\n",(double)ans);
  return 0;
}
