// Independent exhaustive Grassmannian search. No pairing/class formula used.
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
using Row=std::array<int,6>;
int mod(long long x,int p){return (x%p+p)%p;}
int power(int x,int n,int p){int a=1;while(n--){a=a*x%p;}return a;}
int generator(int p){for(int g=2;g<p;++g){int x=1,k=0;do{x=x*g%p;++k;}while(x!=1);if(k==p-1)return g;}return 0;}
int main(int argc,char**argv){
 if(argc!=2)return 2;
 std::ofstream out(argv[1]);
 out<<"p,n0,n1,n2,all_subspaces,pure_survivors,planes,seconds\n";
 for(int p: {7,13,19}){
  int g=generator(p);
  for(int n0=0;n0<=6;++n0)for(int n1=0;n1<=6-n0;++n1){
   int n2=6-n0-n1;
   auto ns=std::array<int,3>{n0,n1,n2};
   if(std::array<int,3>{n1,n2,n0}<ns||std::array<int,3>{n2,n0,n1}<ns)continue;
   if(p==19 && ns!=std::array<int,3>{1,1,4}&&ns!=std::array<int,3>{2,2,2})continue;
   Row coeff{};for(int i=0;i<6;++i)coeff[i]=i<n0?1:(i<n0+n1?g:g*g%p);
   long long all=0,survivors=0,planes=0;
   auto start=std::chrono::steady_clock::now();
   for(int u=0;u<4;++u)for(int v=u+1;v<5;++v)for(int w=v+1;w<6;++w){
    std::array<int,3> piv{u,v,w};std::array<std::vector<Row>,3> opts;
    long long chart=1;
    for(int r=0;r<3;++r){
     std::vector<int> free;
     for(int c=piv[r]+1;c<6;++c)if(c!=u&&c!=v&&c!=w)free.push_back(c);
     int count=1;for(auto c:free)count*=p;chart*=count;
     for(int mask=0;mask<count;++mask){
      Row row{};row[piv[r]]=1;int t=mask;
      for(auto c:free){row[c]=t%p;t/=p;}
      long long pure=0;for(int j=0;j<6;++j)pure+=coeff[j]*power(row[j],3,p);
      if(mod(pure,p)==0)opts[r].push_back(row);
     }
    }
    all+=chart;
    for(auto &a:opts[0])for(auto &b:opts[1])for(auto &c:opts[2]){
     ++survivors;
     std::array<Row,3> rows{a,b,c};bool ok=true;
     for(int i=0;i<3&&ok;++i)for(int k=0;k<3&&ok;++k)if(i!=k){
      long long val=0;
      for(int j=0;j<6;++j)val+=coeff[j]*rows[i][j]*rows[i][j]*rows[k][j];
      ok=mod(val,p)==0;
     }
     if(!ok)continue;
     long long val=0;for(int j=0;j<6;++j)val+=coeff[j]*a[j]*b[j]*c[j];
     if(mod(val,p)==0)++planes;
    }
   }
   double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
   out<<p<<','<<n0<<','<<n1<<','<<n2<<','<<all<<','<<survivors<<','<<planes<<','<<secs<<'\n';out.flush();
   std::cout<<p<<" ("<<n0<<','<<n1<<','<<n2<<") planes="<<planes<<" survivors="<<survivors<<'\n';
  }
 }
}
