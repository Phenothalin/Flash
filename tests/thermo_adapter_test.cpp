// #include "thermo_adapter.hpp"            // 包含 ThermoModel 接口
// #include <iostream>
// #include <cassert>
// #include <vector>

// using namespace thermo;

// /// 一个简单的 FakeThermoModel，用于测试接口连通性
// struct FakeThermoModel : ThermoModel {
//   std::vector<double> chemicalPotentials(const PhaseState& st) override {
//     // 返回每个成分 μ = temperature * 1.0 + pressure * 0.0 + n_i * 0.0
//     return std::vector<double>(st.molefraction.size(), st.temperature);
//   }
//   std::vector<std::vector<double>> dMu_dN(const PhaseState& st) override {
//     size_t C = st.molefraction.size();
//     // 返回一个单位矩阵 * 0.5
//     std::vector<std::vector<double>> J(C, std::vector<double>(C, 0.0));
//     for (size_t i = 0; i < C; ++i) J[i][i] = 0.5;
//     return J;
//   }
//   std::vector<double> dMu_dT(const PhaseState& st) override {
//     return std::vector<double>(st.molefraction.size(), 1.0);
//   }
//   std::vector<double> dMu_dP(const PhaseState& st) override {
//     return std::vector<double>(st.molefraction.size(), 2.0);
//   }
// };

// int main() {
//   // 构造 FakeThermoModel
//   FakeThermoModel model;

//   // 定义测试相态
//   PhaseState state;
//   state.temperature = 300.0;
//   state.pressure    = 101325.0;
//   state.molefraction = {1.0, 2.0, 3.0};

//   // 测试 chemicalPotentials
//   auto mu = model.chemicalPotentials(state);
//   assert(mu.size() == 3);
//   for (double v : mu) {
//     assert(v == 300.0);
//   }
//   std::cout << "chemicalPotentials OK\n";

//   // 测试 dMu_dN
//   auto J = model.dMu_dN(state);
//   assert(J.size() == 3 && J[0].size() == 3);
//   for (size_t i = 0; i < J.size(); ++i) {
//     for (size_t k = 0; k < J.size(); ++k) {
//       double expected = (i==k ? 0.5 : 0.0);
//       assert(J[i][k] == expected);
//     }
//   }
//   std::cout << "dMu_dN OK\n";

//   // 测试温度／压力导数
//   auto dT = model.dMu_dT(state);
//   auto dP = model.dMu_dP(state);
//   assert(dT[0] == 1.0 && dP[0] == 2.0);
//   std::cout << "dMu_dTemperature / dMu_dPressure OK\n";

//   std::cout << "ThermoModel interface test passed!\n";
//   return 0;
// }
