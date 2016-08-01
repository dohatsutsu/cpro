#include<bits/stdc++.h>
using namespace std;
typedef long double ldouble;

typedef complex<ldouble> P;
ldouble cross(P a,P b){ return imag(b*conj(a)); }
ldouble dot(P a,P b){ return real(b*conj(a)); }
ldouble Arg(P a,P b){ return arg(b*conj(a)); }

ldouble eps=1e-8;
ldouble PI=acos(-1.0);

bool eq(ldouble a, ldouble b){
  return (-eps < a-b && a-b < eps );
}

ldouble Sqrt(ldouble x){
  return sqrt( abs(x) );
}
  
struct S{
  P p,q;
  S(P p,P q):p(p),q(q){}
};
  
struct C{
  P p;
  ldouble r;
  C(P p,ldouble r):p(p),r(r){}
  C(ldouble x,ldouble y,ldouble r):p(P(x,y)),r(r){}
};
//28.983094894 28.983094894 28.983094894 32.013426145 61.883234457 25.958852614
//59 99 2639 2117 670

vector<P> intersect(C a,C b){
  vector<P> res;
  P base=b.p-a.p;
  
  ldouble dist=abs(base);
  if(dist>a.r+b.r+eps)return res;
  if(dist<abs(a.r-b.r)-eps)return res;
  
  ldouble ca = (a.r*a.r+dist*dist-b.r*b.r) / (2.0*a.r*dist);
  ldouble sa = Sqrt(1.0-ca*ca);
  
  P fi=a.p+base/dist*a.r*P(ca,sa);
  P se=a.p+base/dist*a.r*P(ca,-sa);
  
  res.push_back(fi);
  res.push_back(se);
  return res;
}


ldouble intersectArea(C a,C b){
  ldouble dist=abs(a.p-b.p);
  if(dist>a.r+b.r+eps)return 0;
  if(dist<abs(a.r-b.r)+eps)return min(a.r*a.r,b.r*b.r)*PI;
  vector<P> v=intersect(a,b);

  if(v.size()<2)return 0;
  
  ldouble sa=0,sb=0;
  sa=a.r*a.r*PI*abs(Arg(v[0]-a.p,v[1]-a.p))/(2.0*PI);
  sa-=abs( cross(v[0]-a.p,v[1]-a.p))/2.0;

  sb=b.r*b.r*PI*abs(Arg(v[0]-b.p,v[1]-b.p))/(2.0*PI);
  sb-=abs( cross(v[0]-b.p,v[1]-b.p))/2.0;
  
  if( dot(b.p-a.p,v[0]-a.p) < 0 ){
    sa=a.r*a.r*PI-sa;
  }
  if( dot(a.p-b.p,v[0]-b.p) < 0 ){
    sb=b.r*b.r*PI-sb;
  }
  
  return sa+sb;
}
  
  
int main(){
  /*
  double px,py,pr;
  cin>>px>>py>>pr;
  C ca=C(px,py,pr);
  cin>>px>>py>>pr;
  C cb=C(px,py,pr);
  cout<< (double)abs(ca.p-cb.p) <<endl;
  cout<<(double)intersectArea(ca,cb) <<endl;
  */
  while(1){
    ldouble H,W;
    ldouble A,B,AB;
    cin>>W>>H>>A>>B>>AB;
    if(H==0&&W==0)break;

    bool flg=false;
    if(A<B){
      flg=true;
      swap(A,B);
    }
    
    ldouble ar= sqrt(A/PI);
    ldouble br= sqrt(B/PI);

    if( max(ar,br) * 2 > min(H,W) + eps ){
      cout<<"impossible"<<endl;
      continue;
    }
    
    ldouble L=abs(ar-br),R=ar+br,M;

    assert( eq(0, intersectArea( C(0,0,ar) , C(R,0,br) ) ) );
    assert( eq( br*br*PI , intersectArea( C(0,0,ar) , C(L,0,br) ) ) );
    for(int i=0;i<100;i++){
      M=(L+R)/2.0;
      ldouble val=intersectArea( C(0,0,ar) , C(M,0,br) );
      if(val <= AB+eps )R=M;
      else L=M;
    }


    P pa=P(ar,ar);
    P pb=P(W-br,H-br);
    if( eq(   AB , intersectArea( C(pa,ar) , C(pb,br) ) ) ){
      if(flg){
        swap(pa,pb);
        swap(ar,br);
      }
        
      printf("%.8f %.8f %.8f %.8f %.8f %.8f\n",
             (double)real(pa),(double)imag(pa),(double)ar,
             (double)real(pb),(double)imag(pb),(double)br);

      continue;
    }
    
    P base=pb-pa;
    base/=abs(base);
    base*=L;
    pb=pa+base;

    if( real(pb)+br < W+eps && imag(pb)+br < H+eps &&
        -eps < real(pb)-br  && imag(pb)-br > -eps ){

      if(flg){
        swap(pa,pb);
        swap(ar,br);
      }
        
      printf("%.8f %.8f %.8f %.8f %.8f %.8f\n",
             (double)real(pa),(double)imag(pa),(double)ar,
             (double)real(pb),(double)imag(pb),(double)br);
    }else{
      cout<<"impossible"<<endl;
    }

    
  }
  return 0;
}
