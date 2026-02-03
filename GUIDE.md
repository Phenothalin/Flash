# 使用指南

本指南帮助你从零开始配置Flash项目的开发环境并成功运行。

## 1. 系统要求

- **操作系统**：macOS (arm64) 或 Windows (x64)
- **编译器**：支持C++17的编译器（Clang 10+, GCC 9+, MSVC 2019+）
- **CMake**：3.18或更高版本
- **Git**：用于克隆仓库

## 2. 安装vcpkg

vcpkg是C++包管理器，用于安装项目依赖。

### macOS

```bash
# 克隆vcpkg
cd ~
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg

# 编译vcpkg
./bootstrap-vcpkg.sh

# 设置环境变量（添加到~/.zshrc或~/.bashrc）
echo 'export VCPKG_ROOT=~/vcpkg' >> ~/.zshrc
source ~/.zshrc
```

### Windows

```powershell
# 克隆vcpkg
cd C:\
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg

# 编译vcpkg
.\bootstrap-vcpkg.bat

# 设置环境变量
[Environment]::SetEnvironmentVariable("VCPKG_ROOT", "C:\vcpkg", "User")
```

## 3. 克隆项目

```bash
git clone <repository-url>
cd Flash
```

## 4. 构建项目

```bash
# 创建构建目录
mkdir -p build && cd build

# 配置（vcpkg会自动安装依赖）
cmake ..

# 编译
make
```

首次构建时，vcpkg会自动安装以下依赖：
- Eigen3（线性代数库）
- GTest（单元测试框架）
- spdlog（日志库）

## 5. 运行测试

```bash
# 从build目录运行

# 两相闪蒸测试
./RAND/tests/2phase_tests

# 三相闪蒸测试
./RAND/tests/3phase_rand_test

# ThermoPack集成测试
./Thermo/tests/tp_thermo_check
```

## 6. 使用示例

### 基本两相闪蒸

```cpp
#include "rand_flash.hpp"
#include "thermopack_adapter.hpp"

int main() {
    // 创建热力学后端（Peng-Robinson状态方程）
    auto backend = std::make_unique<ThermoAdapterTP>("C1,C2,C3", "PR");
    auto linSolver = ls::createEigenSolver();
    RandFlash flash(*backend, *linSolver);

    // 设置输入条件
    FlashInput input;
    input.temperature = 300.0;  // K
    input.pressure = 1e6;       // Pa
    input.feedMoles = {0.5, 0.3, 0.2};

    // 求解
    SolveOptions options;
    options.enable_stability_test = true;
    auto result = flash.solve(input, options);

    // 输出结果
    if (result.success) {
        for (size_t j = 0; j < result.numPhases(); ++j) {
            std::cout << "相 " << j << ": beta = " << result.beta(j) << std::endl;
        }
    }
    return 0;
}
```

### 编译自定义程序

```bash
# 在build目录下，假设你的程序是my_flash.cpp
g++ -std=c++17 \
    -I../Thermo/include \
    -I../RAND/include \
    -I../PhaseStability/include \
    -I../external/thermopack/include \
    my_flash.cpp \
    -L. -lrandflash -lphase_stability -lthermo \
    -L../external/thermopack/lib/macos -lthermopack \
    -o my_flash
```

## 7. 常见问题

### Q: cmake报错找不到VCPKG_ROOT

确保已正确设置环境变量：
```bash
echo $VCPKG_ROOT  # 应输出vcpkg路径
```

### Q: 运行时找不到libthermopack.dylib

设置动态库路径：
```bash
export DYLD_LIBRARY_PATH=/path/to/Flash/external/thermopack/lib/macos:$DYLD_LIBRARY_PATH
```

### Q: vcpkg安装依赖失败

尝试手动安装：
```bash
$VCPKG_ROOT/vcpkg install eigen3 gtest spdlog
```

## 8. 下一步

- 阅读[CLAUDE.md](CLAUDE.md)了解详细的开发文档
- 查看`RAND/tests/`目录下的测试用例学习更多用法
- 参考`RAND/include/rand_flash.hpp`了解完整API
