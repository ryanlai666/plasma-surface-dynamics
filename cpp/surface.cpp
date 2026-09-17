// C++17 numerical backend. ABI: phases contain [duration, a, d, c, s, g].
// Output per phase: [time, coverage, removed, deposited, net, roughness, events].
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
#define API extern "C" __declspec(dllexport)
#else
#define API extern "C"
#endif

namespace {
bool valid(const double* phases, int count, int cycles, double layer, double* out) {
    if (!phases || !out || count < 1 || cycles < 1 || !(layer > 0) || !std::isfinite(layer)) return false;
    for (std::int64_t i = 0; i < std::int64_t(count)*6; ++i)
        if (phases[i] < 0 || !std::isfinite(phases[i])) return false;
    return true;
}

// Constant-time uniform selection and state changes via a permutation partition:
// ids[0:nmodified] are modified, ids[nmodified:sites] are bare.
struct Surface {
    std::vector<int> ids, positions;
    std::vector<std::int64_t> heights;
    int modified = 0;
    explicit Surface(int n) : ids(n), positions(n), heights(n, 0) {
        std::iota(ids.begin(), ids.end(), 0);
        std::iota(positions.begin(), positions.end(), 0);
    }
    void swap_positions(int a, int b) {
        std::swap(ids[a], ids[b]);
        positions[ids[a]] = a;
        positions[ids[b]] = b;
    }
    void set(int site, bool state) {
        int pos = positions[site];
        if (state && pos >= modified) {
            swap_positions(pos, modified++);
        } else if (!state && pos < modified) {
            swap_positions(pos, --modified);
        }
    }
};
}

API int surface_abi_version() { return 1; }

API int surface_mean_field(const double* phases, int count, int cycles, double layer, double* out) {
    if (!valid(phases, count, cycles, layer, out)) return 1;
    double theta = 0, removed = 0, deposited = 0, elapsed = 0;
    std::size_t row = 0;
    for (int cycle = 0; cycle < cycles; ++cycle) {
        for (int j = 0; j < count; ++j) {
            const double* r = phases + std::size_t(j)*6;
            const double t = r[0], k = r[1]+r[2]+r[3]+r[4]+r[5];
            double integral = theta*t;
            if (k > 0) {
                double equilibrium = r[1]/k;
                integral = equilibrium*t + (theta-equilibrium)*(-std::expm1(-k*t))/k;
                theta = equilibrium + (theta-equilibrium)*std::exp(-k*t);
            }
            removed += (r[3]*integral + r[4]*t)*layer;
            deposited += r[5]*t*layer;
            elapsed += t;
            double* o = out + row++*7;
            o[0]=elapsed; o[1]=theta; o[2]=removed; o[3]=deposited;
            o[4]=removed-deposited; o[5]=0; o[6]=0;
        }
    }
    return 0;
}

API int surface_kmc(const double* phases, int count, int cycles, double layer,
                    int sites, std::uint64_t seed, std::uint64_t max_events, double* out) {
    if (!valid(phases, count, cycles, layer, out) || sites < 1 || max_events < 1) return 1;
    try {
        Surface surface(sites);
        std::mt19937_64 rng(seed);
        auto uniform = [&]() { return std::generate_canonical<double, 53>(rng); };
        auto choose = [&](int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(rng); };
        std::uint64_t removed=0, deposited=0, events=0;
        double elapsed=0;
        std::size_t row=0;
        for (int cycle=0; cycle<cycles; ++cycle) {
            for (int j=0; j<count; ++j) {
                const double* r=phases+std::size_t(j)*6;
                double time=0;
                while (time < r[0]) {
                    const int n=surface.modified;
                    double h[5]={r[1]*(sites-n),r[2]*n,r[3]*n,r[4]*sites,r[5]*sites};
                    const double total=std::accumulate(h,h+5,0.0);
                    if (total == 0) break;
                    if (!std::isfinite(total)) return 1;
                    // 1-U is in (0,1], so log never receives zero.
                    time += -std::log1p(-uniform())/total;
                    if (time > r[0]) break;
                    if (++events > max_events) return 2;
                    const double target=uniform()*total;
                    double cumulative=0;
                    int event=4;
                    for (int e=0; e<5; ++e) {
                        cumulative += h[e];
                        if (target < cumulative) { event=e; break; }
                    }
                    int site;
                    if (event == 0) site=surface.ids[choose(n,sites-1)];
                    else if (event < 3) site=surface.ids[choose(0,n-1)];
                    else site=choose(0,sites-1);
                    surface.set(site,event==0);
                    if (event==2 || event==3) { --surface.heights[site]; ++removed; }
                    if (event==4) { ++surface.heights[site]; ++deposited; }
                }
                elapsed += r[0];
                double mean=std::accumulate(surface.heights.begin(),surface.heights.end(),0.0)/sites;
                double variance=0;
                for (auto h : surface.heights) variance+=(h-mean)*(h-mean);
                double* o=out+row++*7;
                o[0]=elapsed; o[1]=double(surface.modified)/sites;
                o[2]=double(removed)*layer/sites; o[3]=double(deposited)*layer/sites;
                o[4]=(double(removed)-double(deposited))*layer/sites;
                o[5]=std::sqrt(variance/sites)*layer; o[6]=double(events);
            }
        }
    } catch (const std::exception&) { return 3; }
    return 0;
}
