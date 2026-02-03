# CLAUDE.md

本文件为Claude Code (claude.ai/code)在处理本仓库代码时提供指导。

## 项目概述

Flash是一个C++17热力学闪蒸计算引擎，用于多相平衡计算。它使用Newton-Raphson迭代和线搜索求解多组分、多相系统的Gibbs自由能最小化问题。

## 构建系统

### 环境准备

构建前需设置以下环境变量：

```bash
export VCPKG_ROOT=/path/to/vcpkg
```

注意：ThermoPack现已内嵌于项目的`external/thermopack/`目录下，不再需要设置`THERMOPACK_DIR`。

### 构建命令

```bash
# 配置和构建
cd /Users/madao/code/Flash
mkdir -p build && cd build
cmake ..
make

# 构建特定目标
make <target_name>

# 清理构建
rm -rf build && mkdir build && cd build && cmake .. && make
```

### 可用构建目标

- `thermo` - 热力学后端库
- `phase_stability` - 相稳定性分析库
- `randflash` - 主闪蒸求解器库
- `compressibility_factor_test` - 状态方程验证测试
- `tp_thermo_check` - ThermoPack集成测试
- `2phase_tests` - 两相闪蒸GTest测试套件
- `2phase_complex_test` - 复杂两相场景测试
- `3phase_rand_test` - 三相闪蒸测试

### 运行测试

```bash
# 从build目录运行
./Thermo/tests/compressibility_factor_test
./Thermo/tests/tp_thermo_check
./RAND/tests/2phase_tests
./RAND/tests/2phase_complex_test
./RAND/tests/3phase_rand_test
```

## 架构设计

### 分层架构

代码库遵循严格的分层架构，依赖关系清晰：

```
RAND (闪蒸求解器)
    ↓ 依赖于
PhaseStability (稳定性分析)
    ↓ 依赖于
Thermo (热力学后端)
    ↓ 依赖于
外部库 (ThermoPack, Eigen3)
```

### 内嵌ThermoPack

ThermoPack内嵌于`external/thermopack/`，结构如下：

```
external/thermopack/
├── CMakeLists.txt          # 创建导入的thermopack目标
├── include/
│   └── cppThermopack/      # C++头文件（11个文件）
│       ├── cubic.h
│       ├── thermo.h
│       └── ...
└── lib/
    └── macos/              # 平台特定库
        ├── libthermopack.dylib
        └── libthermopack.a
```

构建选项：
- `THERMOPACK_USE_STATIC=OFF`（默认）：动态链接
- `THERMOPACK_USE_STATIC=ON`：静态链接

### 核心模块

**Thermo/** - 热力学性质计算
- `thermo_backend.hpp`：抽象`IThermoBackend`接口，定义热力学操作
  - `chemicalPotentials()`, `dmu_dn()`：化学势及其导数
  - `fugacityCoefficients()`, `lnFugacityCoefficients()`：逸度系数计算
  - `vaporPhaseFlag()`, `liquidPhaseFlag()`：相类型标识符
  - `minGibbsPhaseFlag()`：自动选择Gibbs最小根（ThermoPack的`Phase::mingibbs`）
  - `compressibilityFactor()`：Z = PV/(nRT)，用于相类型判定
- `thermopack_adapter.hpp`：`ThermoAdapterTP`封装ThermoPack的立方状态方程（PR, SRK）
- 提供化学势（通过`chemical_potential_tv`，理想+剩余）、逸度系数及其导数
- 化学势计算：TP → TV转换使用`specific_volume`，然后用`PropertyFlag::total`调用`chemical_potential_tv`

**PhaseStability/** - 相稳定性分析
- `phase_stability.hpp/cpp`：`PhaseStabilityAnalyzer`类
- 通过Michelsen逐次替代法实现切平面距离（TPD）最小化
- 多种子点方法检测潜在相
- 支持双参考相（类气相和类液相）

**RAND/** - 多相闪蒸求解器
- `rand_flash.hpp/cpp`：主`RandFlash`类，包含Newton-Raphson求解器
- `rand_init.cpp`：两相、三相和N相系统的初始化例程
- `rand_solver.cpp`：高层求解器入口，带`SolveOptions`控制
- `linear_solver.hpp`：线性代数操作的抽象接口
- `eigen.cpp`：基于Eigen的线性求解器实现

### 关键设计模式

1. **依赖注入**：`RandFlash`接收`IThermoBackend`和`LinearSolverInterface`引用
2. **策略模式**：基于相数的多种初始化策略
3. **接口隔离**：抽象接口解耦实现

## 关键数据结构

### FlashInput
定义闪蒸计算的输入：
- `z`：总体组成（摩尔分数）
- `T`：温度
- `P`：压力
- `phase_labels`：相标识符（如"vapor"、"liquid"）

### MultiFlashResult
包含求解结果：
- `phases`：`PhaseContext`向量，每个包含：
  - `state`：`PhaseState`，含T、P、moleNumbers、phaseFlag
  - `x`：组成（摩尔分数）
  - `mu`：化学势
  - `m`、`M`：Jacobian矩阵
- `success`：收敛状态
- `iterations`：迭代次数
- `mu_infinity_norm`：收敛误差指标
- `pressure`、`temperature`：系统条件
- 便捷方法：`beta(j)`、`n_phase(j)`、`numPhases()`

### SolveOptions
控制求解器行为：
- `max_iterations`：最大Newton迭代次数
- `tolerance`：收敛容差
- `enable_stability_analysis`：切换自动相检测
- `verbose`：调试输出控制

## 算法概述

### 相稳定性分析
1. 使用逐次替代法计算试探组成的TPD
2. 测试多个随机种子点以找到所有潜在相
3. 同时使用类气相和类液相参考相
4. 返回稳定相列表及其组成

### 多相闪蒸求解器
1. 初始化相分率和组成（通过`rand_init.cpp`）
2. 将所有相设置为使用`minGibbsPhaseFlag`（ThermoPack自动选择稳定根）
3. 对约化化学势进行Newton-Raphson迭代：
   - 计算Hessian矩阵和梯度
   - 通过切空间投影修正Hessian正定性
   - 求解线性系统得到Newton步长
   - 应用线搜索和步长调整
4. 基于约化化学势差异检查收敛
5. 收敛后：通过压缩因子确定实际相类型（Z > 0.5 → 气相）
6. 重排相序：气相 → 类油液相 → 富水液相
7. 返回相分率、组成和收敛状态

### 反应体系支持

RAND算法通过元素守恒支持反应平衡：

**元素矩阵(A)**：定义每个组分的元素组成
- A[e][i] = 组分i中元素e的原子数
- 维度：E（元素数）× C（组分数）
- 非反应体系：A = 单位矩阵
- 反应体系：A = 实际化学式矩阵

**守恒约束**：
- 元素守恒：∑_j ∑_i A[e][i] × n_i^(j) = b_e（常数）
- 替代非反应体系中的组分守恒
- 通过元素势自动满足反应平衡

**示例**：烃类体系（C1, C2, C3）
```
       C1   C2   C3
  C  [  1    2    3  ]   (碳原子)
  H  [  4    6    8  ]   (氢原子)
```

## 代码导航提示

### 入口点
- `RandFlash::solve()`位于`RAND/src/rand_solver.cpp` - 主多相闪蒸求解器
- `PhaseStabilityAnalyzer::analyze()`位于`PhaseStability/src/phase_stability.cpp` - 稳定性测试
- `ThermoAdapterTP`位于`Thermo/include/thermopack_adapter.hpp` - 热力学计算

### 理解求解器
1. 从`RAND/include/rand_flash.hpp`开始了解公共API
2. 阅读`RAND/tests/2phase_tests.cpp`查看使用示例
3. 研究`RAND/src/rand_flash.cpp`了解Newton-Raphson核心算法
4. 查看`RAND/src/rand_init.cpp`了解初始化策略

### 添加新功能
- 新热力学模型：实现`IThermoBackend`接口
- 新线性求解器：实现`LinearSolverInterface`接口
- 新初始化策略：添加到`rand_init.cpp`
- 新收敛准则：修改`rand_flash.cpp`中的`checkConvergence()`

### 使用反应体系

**元素矩阵构建器**（`RAND/include/element_matrix_builder.hpp`）：

```cpp
#include "element_matrix_builder.hpp"

// 通过化学式定义组分
std::vector<SpeciesFormula> species = {
    parseFormula("C1", "CH4"),    // 甲烷
    parseFormula("C2", "C2H6"),   // 乙烷
    parseFormula("C3", "C3H8")    // 丙烷
};

// 构建元素矩阵
std::vector<std::string> elementNames;
auto A = buildElementMatrix(species, elementNames);
// 返回：A[0] = [1, 2, 3]（碳），A[1] = [4, 6, 8]（氢）
// elementNames = ["C", "H"]

// 从进料组成计算元素摩尔数
std::vector<double> z = {0.5, 0.3, 0.2};
auto elementMoles = computeElementMoles(A, z);
// 返回：[1.7, 5.4]（碳和氢的摩尔数）
```

**使用元素矩阵运行闪蒸**：

```cpp
// 设置后端和求解器
auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
auto linSolver = ls::createEigenSolver();
RandFlash flash(*backend, *linSolver);

// 准备输入
FlashInput input;
input.temperature = 300.0;  // K
input.pressure = 1e5;       // Pa
input.feedMoles = {0.5, 0.3, 0.2};

// 重要：对于反应体系，必须禁用稳定性测试
// 并手动指定相数
SolveOptions options;
options.enable_stability_test = false;  // 稳定性分析不适用
options.forced_phase_count = 1;         // 手动指定相数
auto result = flash.solve(input, A, options);

// 元素守恒自动强制执行
// 验证：computeElementMoles(A, result.phases[j].state.moleNumbers) = elementMoles
```

**反应体系重要注意事项**：
- 组分摩尔数不守恒；仅元素摩尔数守恒
- 相稳定性分析不适用（组分组成通过反应变化）
- 必须使用`enable_stability_test = false`和`forced_phase_count`
- 某些组分初始可能不存在（通过反应生成）

## 当前开发状态

**分支**：`feature/reaction`

最近更改：
- 将ThermoPack内嵌到项目（`external/thermopack/`），消除外部依赖
- 将`thermo_adapter.hpp`重命名为`thermopack_adapter.hpp`以提高清晰度
- 用ThermoPack原生的`chemical_potential_tv`（理想+剩余）替换组装的化学势
- 通过`SolveOptions::enable_stability_analysis`添加可切换的稳定性分析控制
- 统一多相接口（移除单独的`solveMultiPhase`方法）
- 使用约化化学势改进收敛准则
- 将初始化逻辑分离到`rand_init.cpp`
- 添加`minGibbsPhaseFlag()`以在迭代期间自动选择Gibbs最小根
- 精简`MultiFlashResult` - 移除冗余的`n_phase`/`beta`字段，数据现通过`phases`访问
- 使用压缩因子Z进行收敛后相类型判定
- **新增**：通过元素守恒扩展RAND算法以支持反应体系
- **新增**：添加元素矩阵构建器（`element_matrix_builder.hpp/cpp`）用于化学式解析
- **新增**：实现反应体系初始化函数（`initializeReactiveSinglePhase`、`initializeReactiveMultiPhase`）
- **新增**：添加反应体系综合测试套件（`reactive_flash_test.cpp`）
- **新增**：在整个Newton-Raphson迭代过程中自动强制执行元素守恒

## 依赖项

- **ThermoPack**：热力学性质计算（内嵌于`external/thermopack/`）
- **Eigen3**：线性代数库（通过vcpkg安装）
- **GTest**：单元测试框架（通过vcpkg安装）
- **spdlog**：日志库（通过vcpkg安装）
- **vcpkg**：C++包管理器（必须设置`VCPKG_ROOT`）

## 文件命名约定

- 头文件：`.hpp`扩展名
- 实现文件：`.cpp`扩展名
- 测试文件：`*_test.cpp`或`*_tests.cpp`后缀
- 每个模块都有`include/`、`src/`和`tests/`子目录
