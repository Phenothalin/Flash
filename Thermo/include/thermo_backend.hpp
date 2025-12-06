#pragma once
#include <vector>
#include <memory>
#include "thermo_adapter.hpp"
#include <sstream>

namespace thermo {

// 描述“某一相”的状态：T、P、该相的摩尔数向量 n、以及热力学后端的相标志
struct PhaseState {
  double Temperature;
  double Pressure;
  std::vector<double> moleNumbers;  // 相内摩尔数
  int phaseFlag;          // 比如 ThermoPack: LIQPH / VAPPH
};

// 统一的热力学后端接口
class IThermoBackend {
public:
  virtual ~IThermoBackend() = default;

  virtual std::vector<double> chemicalPotentials(const PhaseState& st) const = 0;     
  virtual void wilsonK(double T, double P, std::vector<double>& K) const = 0;

  // d(μ/RT)/dn: nc x nc
  virtual std::vector<std::vector<double>> dmu_dn(const PhaseState& st) const = 0;

  // 两个辅助接口：给出“气相/液相”对应的 phaseFlag（方便两相、将来多相入口用）
  virtual int vaporPhaseFlag() const = 0;
  virtual int liquidPhaseFlag() const = 0;

  // 获取组分名称列表
  virtual std::vector<std::string> getComponentNames() const = 0;
};

// ThermoPack 具体实现：把 ThermoAdapterTP 包起来，塞给 RandFlash 用
class ThermoPackBackend : public IThermoBackend {
  public:
    explicit ThermoPackBackend(const std::string& components_csv,
                               const std::string& eos   = "PR",
                               const std::string& mixing= "vdW",
                               const std::string& alpha = "Classic",
                               const std::string& ref   = "Default",
                               bool volume_shift = false)
      : tp_(components_csv, eos, mixing, alpha, ref, volume_shift)
    {
      parseComponentNames(components_csv);
    }

    
    std::vector<double> chemicalPotentials(const PhaseState& st) const override
    {
      return tp_.chemicalPotentials(st.Temperature, st.Pressure, st.moleNumbers, st.phaseFlag);
    }
    void wilsonK(double T, double P, std::vector<double>& K) const override {
      tp_.wilsonK(T, P, K);
    }
    std::vector<std::vector<double>> dmu_dn(const PhaseState& st) const override
    {
      return tp_.dmu_dn(st.Temperature, st.Pressure, st.moleNumbers, st.phaseFlag);
    }
  
    int vaporPhaseFlag() const override { return tp_.VAPPH(); }
    int liquidPhaseFlag() const override { return tp_.LIQPH(); }
  
    std::vector<std::string> getComponentNames() const override {
      return comp_names_;
    }
  private:
    ThermoAdapterTP tp_;
    std::vector<std::string> comp_names_;

    void parseComponentNames(const std::string& csv) {
      std::stringstream ss(csv);
      std::string item;
      while (std::getline(ss, item, ',')) {
          // 去除首尾空格
          size_t first = item.find_first_not_of(" \t");
          if (std::string::npos == first) {
              // 全是空格或空字符串
              comp_names_.push_back(item); 
          } else {
              size_t last = item.find_last_not_of(" \t");
              comp_names_.push_back(item.substr(first, (last - first + 1)));
          }
      }
  }
};

} // namespace rf
