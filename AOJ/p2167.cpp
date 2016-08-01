/*
  AOJ 2167 Find the Point

  Nつの直線が与えられるので、どの直線とも距離が等しくなっている点を探す問題。
  答えが一意に定まらないときはManyと出力する。
  存在しない場合はNoneと出力する。
  
  解法
  直線の数が2つ以下のときは、絶対にManyを出力する。
  
  3つ以上あるとき、まず最初の3つ直線に対してのみ距離が等しくなっている点の集合求める。、
  あとはその点1つ1つに対して、残りの直線との距離も等しくなっているかどうかを確認する。
 */



#include<bits/stdc++.h>
using namespace std;
#define MAX_N 105
typedef complex<double> P;

double eps=1e-6;

bool eq(double a,double b){
  
  return (-eps < a-b && a-b < eps);
}

struct S{
  P s,t;
};

P intersect(P a,P b,P c,P d){
  a-=d,b-=d,c-=d;
  return d+a+(b-a)*imag(a/c)/imag(a/c-b/c);
}

P intersect(S a,S b){
  return intersect(a.s,a.t,b.s,b.t);
}

bool isParallel(S a,S b){
  P ap=a.t-a.s;
  P bp=b.t-b.s;
  return eq( 0 , imag( ap/bp ) );
}

double distance(S a,P p){
  return imag( (p-a.s)*conj(a.t-a.s) )/abs(a.t-a.s);
}

int n;
S t[MAX_N];

vector<S> calc(S a,S b){
  vector<S> res;

  P ap=a.t-a.s;
  P bp=b.t-b.s;
  
  if( isParallel(a,b) ){
    P o=(a.s+b.s)*0.5;
    res.push_back( (S){o , o+ap } );
    return res;
  }
  
  P base=intersect(a,b);

  ap/=abs(ap);
  bp/=abs(bp);
  res.push_back( (S){ base,base+ap+bp  });
  res.push_back( (S){ base,base+ap-bp  });
  return res;
}

void solve(){
  if(n<=2){
    cout<<"Many"<<endl;
    return;
  }
  vector< P > vec;
  
  S a=t[0],b=t[1],c=t[2];
  vector< S > va=calc(a,b),vb=calc(b,c),vc=calc(c,a);
  for(int i=0;i<(int)va.size();i++){
    for(int j=0;j<(int)vb.size();j++){
      for(int k=0;k<(int)vc.size();k++){
        S ab=va[i];
        S bc=vb[j];
        S ca=vc[k];
        if( isParallel(ab,bc) || isParallel(bc,ca) || isParallel(ca,ab) ){
          continue;
        }


        P target=intersect(ab,bc);
        P q0=intersect(bc,ca);
        P q1=intersect(ca,ab);
        if( abs(target-q0) > eps )continue;
        if( abs(target-q1) > eps )continue;
        bool flg=true;
        double dist= abs( distance( t[0] , target) );
        for(int id=0;id<n;id++){
          if( !eq( dist, abs(distance( t[id] , target) ) ) )
            flg=false;
        }

        if(flg)vec.push_back(target);
        
      }
    }
  }
  vector< P > ans;
  for(int i=0;i<(int)vec.size();i++){
    bool flg=true;
    for(int j=0;j<(int)ans.size();j++){
      if( abs(vec[i]-ans[j]) < eps )flg=false;
    }
    if(flg)ans.push_back(vec[i]);
  }

  
  if(ans.size()==0){
    cout<<"None"<<endl;
  }else if(ans.size()==1){
    printf("%.10f %.10f\n",real(ans[0]),imag(ans[0]));
  }else{
    cout<<"Many"<<endl;
  }
}

int main(){
  while(1){
    cin>>n;
    if(n==0)break;
    for(int i=0;i<n;i++){
      double x,y;
      cin>>x>>y;
      t[i].s=P(x,y);
      cin>>x>>y;
      t[i].t=P(x,y);


      
    }
    solve();
  }
  return 0;
}
