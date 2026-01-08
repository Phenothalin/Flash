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

int main(){
  // 1) 选择 EoS 与物系（示例：天然气常见组分）
  // 构造函数签名见 cubic.h: Cubic(comps, eos, mixing="vdW", alpha="Classic", ref="Default", volume_shift=false)

  std::string comps ="N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7";
  Cubic eos("C2,C3","PR", "vdW", "Classic", "Default", false);
  Cubic eos2(comps,"PR", "vdW", "Classic", "Default", false);
  std::vector<double> z_test = { 0.1, 0.015, 0.55, 0.14, 0.12,0.05, 0.045, 0.03, 0.025, 0.012, 0.01};
  std::vector<double> z_normalized_sum = {0.0911577, 0.0136737, 0.501367, 0.127621, 0.109389, 0.0455788, 0.0410209, 0.0273473, 0.0227894, 0.0109389, 0.00911577};
  // 2) 条件与组成
  const double T = 295.0;            // K
  const double p = 101325;            // Pa
  std::vector<double> z = {0.5,0.5};
  std::vector<double> x0 = normalize(z);

  auto prop_vap = eos2.thermo(T, p, z_normalized_sum, eos.VAPPH, /*dlnfugdt*/false, /*dlnfugdp*/false, /*dlnfugdn*/false);
  auto prop_liq = eos2.thermo(T, p, z_normalized_sum, eos.LIQPH, /*dlnfugdt*/false, /*dlnfugdp*/false, /*dlnfugdn*/false);
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
  
  auto FlashResult2 = eos2.two_phase_tpflash(295, 2e6, z_normalized_sum);
  auto FlashResult = eos.two_phase_tpflash(295, 101325, z);
  std::cout << "Vapor fraction: " << FlashResult.betaV << std::endl;
  std::cout << "Liquid composition: ";
  for (const double& xi : FlashResult.x) std::cout << xi << " ";
  std::cout << std::endl;
  std::cout << "Vapor composition: ";
  for (const double& yi : FlashResult.y) std::cout << yi << " ";
  std::cout << std::endl;
}
