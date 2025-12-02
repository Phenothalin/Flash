// tests/tp_thermo_check.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <cppThermopack/cubic.h>   // ThermoPack C++ 包装头
#include <cppThermopack/thermo.h>  // Thermo 基类与 Phase flags


static inline std::vector<double> normalize(const std::vector<double>& z){
  double s = 0.0; for(double v : z) s += v;
  std::vector<double> x(z.size());
  for(size_t i=0;i<z.size();++i) x[i] = z[i]/s;
  return x;
}

static inline double max_abs_diff(const std::vector<std::vector<double>>& A,
                                  const std::vector<std::vector<double>>& B){
  double m = 0.0;
  for(size_t i=0;i<A.size();++i)
    for(size_t j=0;j<A[i].size();++j)
      m = std::max(m, std::abs(A[i][j]-B[i][j]));
  return m;
}

// 组装 d(ln x)/dn
static std::vector<std::vector<double>>
make_dlnx_dn(const std::vector<double>& n){
  const size_t nc = n.size();
  double N = 0.0; for(double v:n) N += v;
  std::vector<std::vector<double>> M(nc, std::vector<double>(nc, 0.0));
  for(size_t i=0;i<nc;++i){
    for(size_t j=0;j<nc;++j){
      double delta = (i==j) ? 1.0 : 0.0;
      M[i][j] = (delta/std::max(n[i],1e-300)) - (1.0/N);
    }
  }
  return M;
}

int main(){
  // 1) 选择 EoS 与物系（示例：天然气常见组分）
  // 构造函数签名见 cubic.h: Cubic(comps, eos, mixing="vdW", alpha="Classic", ref="Default", volume_shift=false)
  Cubic eos("N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7",
            "PR", "vdW", "Classic", "Default", false);

  // 2) 条件与组成
  const double T = 295.0;            // K
  const double p = 2.0e6;            // Pa
  std::vector<double> z = {0.03, 0.015, 0.55, 0.14, 0.12,0.05, 0.045, 0.03, 0.025, 0.012, 0.01};
  std::vector<double> x0 = normalize(z);

  auto prop_vap = eos.thermo(T, p, z, eos.VAPPH, /*dlnfugdt*/false, /*dlnfugdp*/false, /*dlnfugdn*/false);
  auto prop_liq = eos.thermo(T, p, z, eos.LIQPH, /*dlnfugdt*/false, /*dlnfugdp*/false, /*dlnfugdn*/false);
  const std::vector<double>& lnphi_vap = prop_vap.value();
  const std::vector<double>& lnphi_liq = prop_liq.value();
  for(size_t i = 0; i < lnphi_vap.size() ; ++i){
    std::cout<<lnphi_vap[i]<<"  ";
  }
  std::cout<<std::endl;
  for(size_t i = 0; i < lnphi_liq.size() ; ++i){
    std::cout<<lnphi_liq[i]<<"  ";
  }
  std::cout<<std::endl;

  auto z_vap = eos.zfac(T, p, z, eos.VAPPH);
  auto z_liq = eos.zfac(T, p, z, eos.LIQPH);
  std::cout << "Z_vap = " << z_vap.value()
            << ", Z_liq = " << z_liq.value() << std::endl;

  // 3) 先做两相 TP-flash，拿到两相组成（以便分别在两相上评估导数）
  auto fr = eos.two_phase_tpflash(T, p, x0);
  auto xL = fr.x;    // 液相组成
  auto yV = fr.y;    // 气相组成

  // 4) 在两相上分别评估 ln(phi) 的导数 dlnphi/dn（Tp 接口）
  //    注意：thermo(...) 的 dlnfugdn 返回的是 ln(逸度系数) 的导数（文档所述）
  //    下面仅取 dn()，不访问任何私有 value_
  auto propL = eos.thermo(T, p, xL, eos.LIQPH, /*dlnfugdt*/false, /*dlnfugdp*/false, /*dlnfugdn*/true);
  auto propV = eos.thermo(T, p, yV, eos.VAPPH, /*dlnfugdt*/false, /*dlnfugdp*/false, /*dlnfugdn*/true);

  const auto& dlnphi_dn_L = propL.dn();  // nc x nc
  const auto& dlnphi_dn_V = propV.dn();

  // 5) 自己补上 d(ln x)/dn，得到 d(ln f)/dn = d(ln phi)/dn + d(ln x)/dn
  //    （p 对 n 不敏感，d(ln p)/dn = 0）
  //    组装时相内 n 可取为 x，因为只差一个总量因子，对 d(ln x)/dn 的表达式只需要 n_i 与 N 的比例。
  std::vector<double> nL = xL, nV = yV;
  auto dlnx_dn_L = make_dlnx_dn(nL);
  auto dlnx_dn_V = make_dlnx_dn(nV);

  // 组装相应的 d(ln f)/dn 矩阵
  const size_t nc = xL.size();
  std::vector<std::vector<double>> dlnf_dn_L(nc, std::vector<double>(nc, 0.0));
  std::vector<std::vector<double>> dlnf_dn_V(nc, std::vector<double>(nc, 0.0));
  for(size_t i=0;i<nc;++i){
    for(size_t j=0;j<nc;++j){
      dlnf_dn_L[i][j] = dlnphi_dn_L[i][j] + dlnx_dn_L[i][j];
      dlnf_dn_V[i][j] = dlnphi_dn_V[i][j] + dlnx_dn_V[i][j];
    }
  }

  // 6) 标准一致性检查：
  //    (a) 伸缩不变性：d(ln f)/dn · n ≈ 0
  auto check_scale = [&](const std::vector<std::vector<double>>& A,
                         const std::vector<double>& n,
                         const char* tag){
    double worst = 0.0;
    for(size_t i=0;i<nc;++i){
      double s = 0.0;
      for(size_t j=0;j<nc;++j) s += A[i][j]*n[j];
      worst = std::max(worst, std::abs(s));
    }
    std::cout << "[scale-invariance] max|row·n| ("<<tag<<") = " << worst << "\n";
  };
  check_scale(dlnf_dn_L, nL, "L");
  check_scale(dlnf_dn_V, nV, "V");

  //    (b) 左右差分对称性（数值扰动检验，轻量）
  auto finite_diff = [&](const std::vector<double>& x, int phase){
    const double eps = 1e-8;
    std::vector<std::vector<double>> J(nc, std::vector<double>(nc, 0.0));
    for(size_t j=0;j<nc;++j){
      std::vector<double> xp=x, xm=x;
      xp[j] = std::max(1e-12, x[j] + eps);
      xm[j] = std::max(1e-12, x[j] - eps);
      // 归一避免漂移
      xp = normalize(xp);
      xm = normalize(xm);
      auto pR = eos.thermo(T, p, xp, phase, false, false, true).dn();
      auto pL = eos.thermo(T, p, xm, phase, false, false, true).dn();
      // 中心差分近似 ∂/∂n_j [ln phi] 的第 i 行：对 dn 再合成 dlnf/dn
      auto dlnx_dn_R = make_dlnx_dn(xp);
      auto dlnx_dn_Lm = make_dlnx_dn(xm);
      for(size_t i=0;i<nc;++i){
        // 取 (lnf)_i 对 n_j 的导数：这里直接用“解析导 + dlnx”，
        J[i][j] = 0.5*((pR[i][j]+dlnx_dn_R[i][j]) + (pL[i][j]+dlnx_dn_Lm[i][j]));
      }
    }
    return J;
  };

  auto Jnum_L = finite_diff(xL, eos.LIQPH);
  auto Jnum_V = finite_diff(yV, eos.VAPPH);

  std::cout << "[consistency] max|dlnf/dn(analytic) - dlnf/dn(numeric)| (L) = "
            << max_abs_diff(dlnf_dn_L, Jnum_L) << "\n";
  std::cout << "[consistency] max|dlnf/dn(analytic) - dlnf/dn(numeric)| (V) = "
            << max_abs_diff(dlnf_dn_V, Jnum_V) << "\n";

  std::cout << "OK\n";
  return 0;
}
