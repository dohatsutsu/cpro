#include<bits/stdc++.h>
using namespace std;

typedef map<int,double> vec;
typedef vector< vec > mat;
double eps=1e-8;
 
int K,N,M;
int x[105],y[105];
int r[3005],s[3005],t[3005];
int C[1005];

//行列のほとんどの要素が0であるため、gaussjordanではなくJacobi法で求めた
vec Jacobi(mat &A,vec &B){
  int N=A.size();
  vec now=B,next;
  
  while(1){
    bool flg=true;
    next.clear();
    for(int i=0;i<N;i++){
      double num=B[i];
      for(vec :: iterator it=A[i].begin();it!=A[i].end();it++)
        if( it->first != i)  num -= (it->second) * now[it->first];

      next[i]=num/A[i][i];

      if( abs(next[i]-now[i]) > (1e-5) )flg=false;
    }
    now=next;
    if(flg)break;
  }
  return now;
}
 
int main(){
  int T,Tc=1;
  scanf("%d",&T);
  while(T--){
    scanf("%d %d %d",&K,&M,&N);
    for(int i=0;i<K;i++){
      scanf("%d %d",&x[i],&y[i]);
    }
 
    for(int i=0;i<M;i++){
      scanf("%d %d %d",&r[i],&s[i],&t[i]);
      if(s[i]>0)s[i]--;
      if(t[i]>0)t[i]--;
    }
 
    mat AX(N),AY(N);
    vec BX,BY;
 
    for(int i=0;i<M;i++){
      int a,b,k;
      a=s[i],b=t[i],k=r[i];
      if(a>b)swap(a,b);
      if(a<0){
        a=(-a)-1;
        AX[b][b]-=k;  AY[b][b]-=k;
        BX[b]-=x[a]*k; BY[b]-=y[a]*k;
      }else{
        AX[a][b]+=k; AY[a][b]+=k;
        AX[b][a]+=k; AY[b][a]+=k;
        AX[b][b]-=k; AY[b][b]-=k;
        AX[a][a]-=k; AY[a][a]-=k;
      }
    }

    
    vec xx=Jacobi(AX,BX);
 
    double ansx[N];
    for(int i=0;i<N;i++){
      ansx[i]=xx[i];
    }
    vec yy=Jacobi(AY,BY);
 
    double ansy[N];
    for(int i=0;i<N;i++){
      ansy[i]=yy[i];
    }
     
    printf("Test case number : %d\n",Tc++);
    for(int i=0;i<N;i++){
      printf("%d %.2f %.2f\n",i+1,ansx[i],ansy[i]);
    }
 
     
  }
  return 0;
}

