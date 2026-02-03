# Flash

高性能C++17热力学闪蒸计算引擎，用于多相平衡计算。

## 项目简介

Flash通过Newton-Raphson迭代和线搜索求解多组分、多相系统的Gibbs自由能最小化问题。支持非反应体系和反应体系的相平衡计算。

## 主要特性

- **多相闪蒸计算**：支持两相、三相和N相平衡
- **相稳定性分析**：基于切平面距离（TPD）最小化的多种子点方法
- **反应体系支持**：通过元素守恒处理化学反应平衡
- **热力学模型**：通过ThermoPack支持Peng-Robinson (PR)、Soave-Redlich-Kwong (SRK)等多种模型
- **鲁棒求解器**：Newton-Raphson迭代配合线搜索和自动相检测
- **现代C++**：清晰的C++17代码库，采用依赖注入和接口隔离设计

## 快速开始

### 环境准备

```bash
# 设置vcpkg根目录
export VCPKG_ROOT=/path/to/vcpkg

# vcpkg会通过vcpkg.json自动安装以下依赖：
# - Eigen3 (线性代数库)
# - GTest (单元测试框架)
# - spdlog (日志库)
```

### 构建

```bash
cd /Users/madao/code/Flash
mkdir -p build && cd build
cmake ..
make
```

### 运行测试

```bash
# 两相闪蒸测试
./RAND/tests/2phase_tests

# 三相闪蒸测试
./RAND/tests/3phase_rand_test

# 反应体系测试
./RAND/tests/reactive_flash_test

# ThermoPack集成测试
./Thermo/tests/tp_thermo_check
```

## 使用示例

### 非反应体系闪蒸

```cpp
#include "rand_flash.hpp"
#include "thermopack_adapter.hpp"

// 设置热力学后端（Peng-Robinson状态方程）
auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
auto linSolver = ls::createEigenSolver();
RandFlash flash(*backend, *linSolver);

// 准备输入
FlashInput input;
input.temperature = 300.0;  // K
input.pressure = 1e6;       // Pa
input.feedMoles = {0.5, 0.3, 0.2};  // C1, C2, C3

// 求解（自动相检测）
SolveOptions options;
options.enable_stability_test = true;
auto result = flash.solve(input, options);

// 访问结果
for (size_t j = 0; j < result.numPhases(); ++j) {
    std::cout << "相 " << j << ": beta = " << result.beta(j) << std::endl;
}
```

### 反应体系闪蒸

```cpp
#include "element_matrix_builder.hpp"

// 通过化学式定义组分
std::vector<SpeciesFormula> species = {
    parseFormula("C1", "CH4"),
    parseFormula("C2", "C2H6"),
    parseFormula("C3", "C3H8")
};

// 构建元素矩阵
std::vector<std::string> elementNames;
auto A = buildElementMatrix(species, elementNames);

// 使用元素守恒求解
SolveOptions options;
options.enable_stability_test = false;  // 反应体系不适用相稳定性分析
options.forced_phase_count = 1;
auto result = flash.solve(input, A, options);
```

## 项目结构

```
Flash/
├── Thermo/              # 热力学后端层
│   ├── include/
│   │   ├── thermo_backend.hpp      # IThermoBackend接口
│   │   └── thermopack_adapter.hpp  # ThermoPack适配器
│   └── tests/
├── PhaseStability/      # 相稳定性分析
│   ├── include/
│   │   └── phase_stability.hpp
│   └── src/
├── RAND/                # 多相闪蒸求解器
│   ├── include/
│   │   ├── rand_flash.hpp          # 主求解器API
│   │   ├── element_matrix_builder.hpp
│   │   └── rand_logger.hpp
│   ├── src/
│   │   ├── rand_flash.cpp          # Newton-Raphson核心
│   │   ├── rand_init.cpp           # 初始化策略
│   │   └── rand_solver.cpp         # 高层求解器
│   └── tests/
└── external/
    └── thermopack/      # 内嵌ThermoPack库
```

## 架构设计

代码库采用分层架构：

```
RAND (闪蒸求解器)
    ↓
PhaseStability (稳定性分析)
    ↓
Thermo (热力学后端)
    ↓
外部库 (ThermoPack, Eigen3)
```

## 算法流程

1. **初始化**：根据相数初始化相分率和组成
2. **相稳定性分析**（可选）：通过TPD最小化检测潜在相
3. **Newton-Raphson迭代**：
   - 从化学势计算Hessian矩阵和梯度
   - 通过切空间投影修正Hessian正定性
   - 求解线性系统得到Newton步长
   - 应用线搜索和步长调整
4. **收敛检查**：基于约化化学势差异
5. **相类型判定**：通过压缩因子判断相类型（Z > 0.5 → 气相）

## 开发指南

### 构建目标

```bash
make thermo                    # 热力学后端库
make phase_stability           # 相稳定性库
make randflash                 # 主闪蒸求解器库
make 2phase_tests              # 两相测试套件
make 3phase_rand_test          # 三相测试
make reactive_flash_test       # 反应体系测试
```

### 添加新功能

- **新热力学模型**：实现`IThermoBackend`接口
- **新线性求解器**：实现`LinearSolverInterface`接口
- **新初始化策略**：添加到`rand_init.cpp`

## 依赖项

- **ThermoPack**：内嵌于`external/thermopack/`（无需外部安装）
- **Eigen3**：线性代数库（通过vcpkg安装）
- **GTest**：单元测试框架（通过vcpkg安装）
- **spdlog**：日志库（通过vcpkg安装）

## 文档

详细开发文档请参见[CLAUDE.md](CLAUDE.md)。

## 许可证

[在此添加许可证信息]
