#include<bits/stdc++.h>
using namespace std;
typedef complex<double> P;
typedef pair< P , P > Pair;
P inf=P( 1e30 , 1e30 );

P reflect(P a,P b,P c){
  b-=a;c-=a;
  return a+conj(c/b)*b;
}

P intersect(P a,P b,P c,P d){
  a-=d;b-=d;c-=d;
  return d+a+(b-a)*imag(-a/c)/imag(b/c-a/c);
}

string s;
int p,len;

Pair calc(Pair a,Pair b){
  if(a.second==inf && b.second==inf){
    return Pair( a.first , b.first );
  }
  if(a.second==inf){
    return Pair( reflect( b.first , b.second , a.first ) , inf );
  }
  if(b.second==inf){
    return Pair( reflect( a.first , a.second , b.first ) , inf );
  }
  return Pair( intersect( a.first , a.second , b.first , b.second ) , inf );
}

Pair solve();

Pair func(){
  if(s[p]=='(' && s[p+1]=='('){
    p++;
    Pair res=solve();
    p++;
    return res;
  }else{
    p++;
    int x=0,y=0;

    bool xf=false;
    if(s[p]=='-')xf=true,p++;
    while( '0'<=s[p] && s[p]<='9'){
      x=x*10+(s[p]-'0');
      p++;
    }
    if(xf)x*=-1;
    
    p++;

    bool yf=false;
    if(s[p]=='-')yf=true,p++;
    while( '0'<=s[p] && s[p]<='9'){
      y=y*10+(s[p]-'0');
      p++;
    }
    if(yf)y*=-1;
    p++;
    return Pair( P(x,y) , inf );
  }
}


Pair solve(){
  Pair a=func();
  while(p<len && s[p]=='@' ){
    p++;
    Pair b=func();
    a=calc(a,b);
  }
  return a;
}

int main(){
  while(1){
    cin>>s;
    if(s=="#")break;
    p=0;
    len=s.size();
    Pair pair=solve();
    P ans=pair.first;
    printf("%.8f %.8f\n",real(ans),imag(ans));
  }
  return 0;
}
