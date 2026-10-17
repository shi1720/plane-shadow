#include <array>
#include <vector>
#include <fstream>
using Row=std::array<int,6>;
int main(int argc,char**argv){if(argc!=2)return 2;std::ofstream out(argv[1]);out<<"p,pattern,lines\n";const int p=7;
for(int which=0;which<2;++which){
 Row coeff=which?Row{1,1,3,3,2,2}:Row{1,1,1,1,3,2};long long count=0;
 for(int u=0;u<5;++u)for(int v=u+1;v<6;++v){
  std::array<int,2> piv{u,v};std::array<std::vector<Row>,2> opts;
  for(int r=0;r<2;++r){std::vector<int> free;for(int j=piv[r]+1;j<6;++j)if(j!=v)free.push_back(j);
   int num=1;for(auto j:free)num*=p;
   for(int mask=0;mask<num;++mask){Row a{};a[piv[r]]=1;int t=mask;
    for(auto j:free){a[j]=t%p;t/=p;}int s=0;for(int j=0;j<6;++j)s+=coeff[j]*a[j]*a[j]*a[j];
    if(s%p==0)opts[r].push_back(a);
   }
  }
  for(auto &a:opts[0])for(auto &b:opts[1]){int s=0,t=0;for(int j=0;j<6;++j){s+=coeff[j]*a[j]*a[j]*b[j];t+=coeff[j]*a[j]*b[j]*b[j];}if(s%p==0&&t%p==0)++count;}
 }
 out<<7<<','<<(which?'Y':'X')<<','<<count<<'\n';
}}
