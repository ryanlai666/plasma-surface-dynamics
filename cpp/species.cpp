// Generic count-based SSA for the explicitly parameterized species network.
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>
#ifdef _WIN32
#define API extern "C" __declspec(dllexport)
#else
#define API extern "C"
#endif
API int species_abi_version(){return 1;}
API int species_kmc(int ns,int ne,int ng,int nt,const std::int64_t* initial,
 const std::int32_t* src,const std::int32_t* dst,const std::int32_t* ads,
 const double* rates,const std::int64_t* delta,const double* times,double exposure,
 std::uint64_t seed,std::uint64_t max_events,std::int64_t* output,
 std::int64_t* gas_output,std::int64_t* event_counts){
 try {
  if(ns<1||ne<1||ng<1||nt<1||!std::isfinite(exposure)||exposure<0)return 1;
  for(int i=0;i<ns;++i)if(initial[i]<0)return 1;
  for(int i=0;i<nt;++i)if(!std::isfinite(times[i])||times[i]<0||(i&&times[i]<times[i-1]))return 1;
  for(int e=0;e<ne;++e)if(src[e]<0||src[e]>=ns||dst[e]<0||dst[e]>=ns||rates[e]<0||!std::isfinite(rates[e]))return 1;
  std::vector<std::int64_t> count(initial,initial+ns),gas(ng,0);std::vector<double> hazard(ne);
  std::mt19937_64 rng(seed);std::uniform_real_distribution<double> uniform(0.,1.);
  double t=0.;int sample=0;std::uint64_t events=0;
  while(sample<nt){
   double total=0.;for(int e=0;e<ne;++e){hazard[e]=((t>=exposure&&ads[e])?0.:rates[e])*count[src[e]];total+=hazard[e];}
   double u=uniform(rng);if(u==0.)u=std::numeric_limits<double>::min();
   double candidate=total>0?t-std::log(u)/total:std::numeric_limits<double>::infinity();
   double boundary=t<exposure?exposure:std::numeric_limits<double>::infinity();double next=candidate<boundary?candidate:boundary;
   while(sample<nt&&times[sample]<=next){for(int i=0;i<ns;++i)output[sample*ns+i]=count[i];for(int g=0;g<ng;++g)gas_output[sample*ng+g]=gas[g];++sample;}
   if(sample==nt)break;
   if(boundary<=candidate){t=boundary;continue;}
   if(++events>max_events)return 2;
   double target=uniform(rng)*total,cumulative=0.;int chosen=-1;
   for(int e=0;e<ne;++e){cumulative+=hazard[e];if(target<cumulative){chosen=e;break;}}
   if(chosen<0||count[src[chosen]]<=0)return 3;
   --count[src[chosen]];++count[dst[chosen]];++event_counts[chosen];
   for(int g=0;g<ng;++g)gas[g]+=delta[chosen*ng+g];
   t=candidate;
  }
  return 0;
 }catch(...){return 3;}
}
