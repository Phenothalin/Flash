// tests/tp_thermo_check.cpp
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <stdexcept>

#include <cppThermopack/thermo.h>
#include <cppThermopack/cubic.h>

static constexpr double Rgas = 8.314462618;

static std::vector<double> normalize(const std::vector<double>& n){
  const double N = std::accumulate(n.begin(), n.end(), 0.0);
  if(N <= 0) throw std::runtime_error("Total moles <= 0");
  std::vector<double> z(n.size());
  for(size_t i=0;i<n.size();++i) z[i] = n[i]/N;
  return z;
}

int main(){
  // 1) 正确创建 PR 模型（构造器里完成 EOS 初始化）
  //    见 cubic.h: PengRobinson(comps, mixing, alpha, ref, volume_shift)
  //    例如混合：N2, CO2, C1, C2, C3, iC4, nC4, iC5, nC5, nC6, nC7
  std::string comps = "N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7";
  PengRobinson model(comps, "vdW", "Classic", "Default", /*volume_shift*/false); // :contentReference[oaicite:2]{index=2}

  const int VAPPH = model.VAPPH; // 相标志在 Thermo 的公有成员里  :contentReference[oaicite:3]{index=3}

  // 2) 状态与物料
  const double T = 295.0;   // K
  const double P = 2.0e6;   // Pa
  std::vector<double> n = {0.000001,0.015,0.55, 0.14,  0.12,  0.05, 0.045, 0.03, 0.025, 0.012, 0.01};
  auto z = normalize(n);
  const size_t C = z.size();

  // 3) 解析导数 d(lnφ)/dn：用 TP 接口 thermo(..., dlnfugdn=true)
  auto dprop = model.thermo(T, P, z, VAPPH,
                            /*dlnfugdt*/false,
                            /*dlnfugdp*/false,
                            /*dlnfugdn*/true);                           // :contentReference[oaicite:4]{index=4}
  const auto& dlnf_dn = dprop.dn(); // vector2d: [i][k]

  // 3.1 尺寸检查
  if(dlnf_dn.size() != C || dlnf_dn[0].size() != C){
    std::cerr << "[tp_thermo_check] ERROR: dlnphi/dn shape mismatch: "
              << dlnf_dn.size() << "x"
              << (dlnf_dn.empty()?0:dlnf_dn[0].size()) << " vs " << C << "x" << C << "\n";
    return 2;
  }

  // 3.2 物理一致性（整体缩放 n 不改变 z）：对每个 i，sum_k d(lnφ_i)/dn_k * n_k ≈ 0
  double max_abs_scale_violation = 0.0;
  for(size_t i=0;i<C;++i){
    double dot = 0.0;
    for(size_t k=0;k<C;++k) dot += dlnf_dn[i][k] * n[k];
    max_abs_scale_violation = std::max(max_abs_scale_violation, std::abs(dot));
  }
  std::cout << "[tp_thermo_check] max |(dlnphi/dn * n)| = "
            << max_abs_scale_violation << " (should be ~ 0)\n";

  // 4) 构造 Rand-Flash 需要的 ∂μ/∂n（不依赖 lnφ 的数值）
  //    μ_i = RT( lnφ_i + ln x_i ) + 常数；因此
  //    ∂μ_i/∂n_k = RT( ∂lnφ_i/∂n_k + ∂ln x_i/∂n_k )
  const double RT = Rgas * T;
  double Ntot = std::accumulate(n.begin(), n.end(), 0.0);
  std::vector<std::vector<double>> dmu_dn(C, std::vector<double>(C, 0.0));
  for(size_t i=0;i<C;++i){
    for(size_t k=0;k<C;++k){
      const double dlnx = ((i==k) ? 1.0/std::max(n[i],1e-300) : 0.0) - 1.0/Ntot;
      dmu_dn[i][k] = RT * ( dlnf_dn[i][k] + dlnx );
    }
  }

  std::cout << "J = dmu/dn (top-left 3x3):\n";
  for(int i=0;i<std::min<size_t>(3,C);++i){
    for(int k=0;k<std::min<size_t>(3,C);++k){
      std::cout << dmu_dn[i][k] << (k+1<(int)std::min<size_t>(3,C) ? ' ' : '\n');
    }
  }
  return 0;
}
