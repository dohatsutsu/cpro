#include<bits/stdc++.h>
using namespace std;
typedef complex<double> P;
typedef vector<P> vecP;
typedef pair<P,P> L;
typedef pair<P,P> S;
typedef pair<P,double> C;
const double eps=1e-8;
const double PI=acos(-1);
const double PI2=PI*2.0;

namespace std{
  bool operator < (const P &a,const P &b){
    return (a.imag()==b.imag()?
            a.real()<b.real():
            a.imag()<b.imag());
  }
};

double Sqrt( double x ){
  if(x<0)return 0;
  else return sqrt(x);
}

P Vector(L a){
  return a.second-a.first;
}

bool eq(double a,double b){
  return (-eps<a-b&&a-b<eps);
}

bool eq(P a,P b){
  return ( eq(a.real(),b.real()) && eq(a.imag(),b.imag()) );
}

double dot(P a,P b){
  return real(b*conj(a));
}

double cross(P a,P b){
  return imag(b*conj(a));
}

double getArg(P a,P b){
  return arg(b*conj(a));
}

P project(P a,P b,P c){
  b-=a,c-=a;
  return a+b*real(c/b);
}

P reflect(P a,P b,P c){
  b-=a,c-=a;
  return a+b*real(c/b);
}

int ccw(P a,P b,P c){
  P ab=b-a,ac=c-a;
  P k=ac*conj(ab);
  if(k.imag()>0)return 1;
  if(k.imag()<0)return -1;
  if(k.real()<0)return 2;
  if(abs(ab)<abs(ac))return -2;
  return 0;
}

bool isParallel(P a,P b){
  return eq(0, cross(a,b));
}

bool isParallel(L l0,L l1){
  return eq(0, cross( Vector(l0) , Vector(l1) ) );
}

bool onLP(L l,P p){
  P a=l.first, b=l.second;
  return eq(0, cross(b-a,p-a));
}

bool onSP(S s,P p){
  P a=s.first, b=s.second;
  return eq( abs(b-a) , abs(a-p)+abs(b-p) );
}

bool isCrossSS(S s0,S s1){
  P a=s0.first, b=s0.second;
  P c=s1.first, d=s1.second;
  int f0 = ccw(a,b,c) * ccw(a,b,d);
  int f1 = ccw(c,d,a) * ccw(c,d,b);
  return (f0<=0 && f1<=0);
}

bool isCrossLS(L l,S s){
  P a=l.first, b=l.second;
  P c=s.first, d=s.second;
  return ( ccw(a,b,c) * ccw(a,b,d) <= 0 );
}

double distLP(L l,P p){
  P a=l.first, b=l.second;
  double res = cross(b-a,p-a) / abs(b-a);
  return abs(res);
}

double distSP(S s,P p){
  P a=s.first, b=s.second;
  if( dot(b-a,p-a) < eps )return abs(p-a);
  if( dot(a-b,p-b) < eps )return abs(p-b);
  return distLP(s,p);
}

P getCrossLL(L l0,L l1){
  P a=l0.first, b=l0.second;
  P c=l1.first, d=l1.second;
  a-=d;b-=d;c-=d;
  return d+a+(b-a)*imag(a/c)/imag(a/c-b/c);
}


 
int inPolygon(vecP &t,P p){
  int n=t.size();
  double sum=0;
  for(int i=0;i<n;i++){
    P a=t[i],b=t[(i+1==n?0:i+1)];
    if( onSP(S(a,b),p) )return 1;
    sum+= getArg(a-p,b-p);
  }
  if( abs(sum) < eps )return 0;
  else return 2;
}

vecP andrewScan(vecP &t){
  int N=t.size(),C=0;
  vecP R(N);
  for(int i=0;i<N;i++){
    while(2<=C&&ccw(R[C-2],R[C-1],t[i])==-1)C--;
    R[C++]=t[i];
  }
  vecP res(C);
  for(int i=0;i<C;i++)res[i]=R[i];
  return res;
}
 
vecP ConvexHull(vecP &t){
  sort(t.begin(),t.end());
  vecP u=andrewScan(t);
  reverse(t.begin(),t.end());
  vecP l=andrewScan(t);
  for(int i=1;i+1<(int)l.size();i++)u.push_back(l[i]);
  return u;
}

vecP cutConvex(vecP &t,L l){
  P a=l.first, b=l.second;
  int N=t.size();
  vecP res;
  for(int i=0;i<N;i++){
    P c=t[i],d=t[(i+1)%N];
    int C=ccw(a,b,c),D=ccw(a,b,d);
    if(C!=-1)res.push_back(c);
    if(C==-D&&abs(C)==1)res.push_back(getCrossLL( l ,L(c,d) ));
  }
  return res;
}

double maxDist(vecP &t){
  vecP u=t;//ConvexHull(t);
  int N=u.size(),K=0;
  double res=0;
  for(int i=0;i<N;i++){
    while(abs(u[i]-u[K])<abs(u[i]-u[(K+1)%N]))K=(K+1)%N;
    res=max(res,abs(u[i]-u[K]));
  }
  return res;
}


bool compare_y(const P &a,const P &b){
  return a.imag() < b.imag();
}

double closest_pair(P *a, int n){
  if(n <= 1) return 1e30;
  int m = n / 2;
  double x = a[m].real();
  double d = min(closest_pair(a, m), closest_pair(a + m, n - m));
  inplace_merge(a, a + m, a + n, compare_y);
  vector<P> b;
  for(int i=0;i<n;i++){
    if( abs(a[i].real() - x) >= d) continue;
    for(int j=0;j<(int)b.size();j++){
      double dx = real(a[i] - b[b.size() - j - 1]);
      double dy = imag(a[i] - b[b.size() - j - 1]);
      if(dy >= d) break;
      d = min(d, sqrt(dx * dx + dy * dy));
    }
    b.push_back(a[i]);
  }
  return d;
}

P _pool[200005];
double minDist(vecP &t){
  int n=t.size();
  for(int i=0;i<n;i++)_pool[i]=t[i];
  sort( _pool, _pool+n);
  return closest_pair(_pool, n);
}

int getStateCC(C a,C b){
  double ar=a.second, br=b.second;
  double dist=abs(a.first-b.first);
  if(dist>ar+br+eps)return 4;
  if(dist>ar+br-eps)return 3;
  if(dist>abs(ar-br)+eps)return 2;
  if(dist>abs(ar-br)-eps)return 1;
  return 0;
}

P getCrossCC(C a,C b){
  P p1=a.first, p2=a.second;
  double r1=a.second, r2=b.second;
  double cA = (r1*r1+norm(p1-p2)-r2*r2) / (2.0*r1*abs(p1-p2));
  return p1+(p2-p1)/abs(p1-p2)*r1*P(cA,Sqrt(1.0-cA*cA));
}

S getTangentCP(C a,P p){
  P base=a.first-p;
  double ar=a.second;
  double w=Sqrt(norm(base)-ar*ar);
  P s=p+base*P(w,ar)/norm(base)*w;
  P t=p+base*P(w,-ar)/norm(base)*w;
  return S(s,t);
}

S getInTangent(C a,C b,double flg=1.0){
  P ap=a.first,bp=b.first;
  double ar=a.second,br=b.second;
  
  P base=bp-ap;
  double w=ar+br;
  double h=Sqrt(norm(base)-w*w);
  P k=base*P(w,h*flg)/norm(base);
  return S(ap+k*ar,bp-k*br);
}
  
S getOutTangent(C a,C b,double flg=1.0){
  P ap=a.first,bp=b.first;
  double ar=a.second,br=b.second;
  
  P base=bp-ap;
  double h=br-ar;
  
  double w=Sqrt(norm(base)-h*h);
  P k=base*P(w,h*flg)/norm(base)*P(0,flg);
  return S(ap+k*ar,bp+k*br);
}
  
vector<S> getTangent(C a,C b){
  P ap=a.first,bp=b.first;
  double ar=a.second,br=b.second;
  
  vector<S> res;
  double dist=abs(ap-bp);
    
  if(dist>ar+br+eps)
    res.push_back(getInTangent(a,b,1));
  
  if(dist>ar+br-eps)
    res.push_back(getInTangent(a,b,-1));
  
  if(dist>abs(ar-br)+eps)
    res.push_back(getOutTangent(a,b,1));
  
  if(dist>abs(ar-br)-eps)
    res.push_back(getOutTangent(a,b,-1));
  
  return res;
}

double getTime(P a,P b){
  assert( eq(cross(a,b),0) );
  return ( dot(a,b) < 0 ? -1.0 : 1.0 ) * abs(b);
}

vecP getCrossCS(C cir,S s, bool debug=false){
  P a=s.first, b=s.second;
  double cr=cir.second;
  P cp=cir.first;
  
  vecP res;
  P base=b-a,  target=project(a,b,cp);
  
  double length=abs(base), h=abs(cp-target);
  base/=length;
  
  if(cr+eps<h)return res;
  double w=Sqrt(cr*cr-h*h);
  double L=getTime(b-a,target-a)-w,  R=L+w*2.0;
  
  if( -eps<L && L< length+eps )res.push_back(a+base*L);
  if( eq(L,R) )return res;
  if( -eps<R && R< length+eps )res.push_back(a+base*R);
  return res;
}
 
double getArea(C c,P a,P b){
  P cp=c.first;
  double cr=c.second;
  
  P va=cp-a,  vb=cp-b;
  double A=abs(va), B=abs(vb);
  double f=cross(va,vb), d=distSP( S(a,b) ,cp), res=0;
  
  if( eq(0, f ) )return 0;
  if(A<cr+eps&&B<cr+eps)return f*0.5;
  if(d>cr-eps)return cr*cr*PI*getArg(va,vb)/PI2;
   
  vecP u=getCrossCS(c, S(a,b) );
  
  assert( !u.empty() );
  u.insert(u.begin(), a),  u.push_back(b);
 
  for(int i=0;i+1<(int)u.size();i++) res+=getArea(c,u[i],u[i+1]);
  return res;
}
 
double getCrossArea(vecP t,C c){
  int n=t.size();
  if(n<3)return 0;
  double res=0;
  for(int i=0;i<n;i++){
    P a=t[i], b=t[(i+1)%n];
    res+=getArea(c,a,b);
  }
  return res;
}

int main(){
  return 0;
}
