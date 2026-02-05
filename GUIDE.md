# Flash 项目 Windows 配置指南

本指南帮助**零基础用户**在全新 Windows 电脑上配置 Flash 项目开发环境。

---

## 目录

1. [环境概述](#1-环境概述)
2. [安装 Git](#2-安装-git)
3. [安装 Visual Studio Build Tools](#3-安装-visual-studio-build-tools)
4. [安装 CMake](#4-安装-cmake)
5. [安装 Ninja](#5-安装-ninja)
6. [安装 vcpkg](#6-安装-vcpkg)
7. [配置环境变量](#7-配置环境变量)
8. [克隆项目](#8-克隆项目)
9. [构建项目](#9-构建项目)
10. [运行测试](#10-运行测试)
11. [常见问题](#11-常见问题)

---

## 1. 环境概述

### 需要安装的工具

| 工具 | 用途 | 版本要求 |
|------|------|----------|
| Git | 版本控制、下载代码 | 任意版本 |
| Visual Studio Build Tools | C++ 编译器 | 2019 或更高 |
| CMake | 构建系统生成器 | 3.18 或更高 |
| Ninja | 快速构建工具 | 任意版本 |
| vcpkg | C++ 包管理器 | 最新版本 |

### 项目依赖库（vcpkg 自动安装）

- **Eigen3** - 线性代数库
- **GTest** - 单元测试框架
- **spdlog** - 日志库

---

## 2. 安装 Git

### 步骤 2.1：下载 Git

1. 打开浏览器，访问：https://git-scm.com/download/win
2. 页面会自动开始下载，如果没有，点击 **"Click here to download manually"**

### 步骤 2.2：安装 Git

1. 双击下载的安装程序（如 `Git-2.xx.x-64-bit.exe`）
2. **所有选项保持默认**，一直点击 **Next**
3. 点击 **Install**，等待安装完成
4. 点击 **Finish**

### 步骤 2.3：验证安装

1. 按 `Win + R`，输入 `cmd`，按回车
2. 输入：
   ```cmd
   git --version
   ```
3. 显示 `git version 2.xx.x` 即成功

---

## 3. 安装 Visual Studio Build Tools

### 步骤 3.1：下载

1. 访问：https://visualstudio.microsoft.com/zh-hans/visual-cpp-build-tools/
2. 点击 **"下载生成工具"**

### 步骤 3.2：安装

1. 运行下载的 `vs_BuildTools.exe`
2. 等待 Visual Studio Installer 加载
3. 勾选 **"使用 C++ 的桌面开发"**
4. 确保右侧选中：
   - MSVC v143 - VS 2022 C++ x64/x86 生成工具
   - Windows 11 SDK（或 Windows 10 SDK）
5. 点击 **"安装"**（需要 10-30 分钟）

### 步骤 3.3：验证

1. 按 `Win` 键，搜索 **"x64 Native Tools Command Prompt"**
2. 打开后输入：
   ```cmd
   cl
   ```
3. 显示 Microsoft C/C++ 编译器版本信息即成功

---

## 4. 安装 CMake

### 步骤 4.1：下载

1. 访问：https://cmake.org/download/
2. 找到 **Windows x64 Installer**
3. 下载 `cmake-3.xx.x-windows-x86_64.msi`

### 步骤 4.2：安装

1. 双击运行 `.msi` 文件
2. 接受许可协议
3. **重要**：选择 **"Add CMake to the system PATH for all users"**
4. 点击 **Install**

### 步骤 4.3：验证

1. **重新打开**命令提示符
2. 输入：
   ```cmd
   cmake --version
   ```
3. 显示 `cmake version 3.xx.x` 即成功

---

## 5. 安装 Ninja

### 步骤 5.1：下载

1. 访问：https://github.com/ninja-build/ninja/releases
2. 下载最新版 `ninja-win.zip`

### 步骤 5.2：安装

1. 解压得到 `ninja.exe`
2. 创建文件夹：`C:\Tools\ninja`
3. 将 `ninja.exe` 复制到该文件夹

> **注意**：Ninja 路径将在第 7 节统一添加到环境变量

---

## 6. 安装 vcpkg

vcpkg 是微软开发的 C++ 包管理器，用于自动安装项目依赖。

### 步骤 6.1：选择安装位置

建议安装在 `C:\vcpkg`（路径不要有中文或空格）

### 步骤 6.2：克隆 vcpkg

1. 打开命令提示符（`Win + R`，输入 `cmd`）
2. 执行以下命令：
   ```cmd
   cd C:\
   git clone https://github.com/microsoft/vcpkg.git
   ```
3. 等待下载完成

### 步骤 6.3：编译 vcpkg

```cmd
cd C:\vcpkg
.\bootstrap-vcpkg.bat
```

等待编译完成，会生成 `vcpkg.exe`

### 步骤 6.4：验证

```cmd
.\vcpkg --version
```

显示版本号即成功

---

## 7. 配置环境变量

需要配置两个环境变量：`VCPKG_ROOT` 和 `PATH`

### 步骤 7.1：打开环境变量设置

1. 按 `Win` 键，搜索 **"环境变量"**
2. 点击 **"编辑系统环境变量"**
3. 点击 **"环境变量"** 按钮

### 步骤 7.2：添加 VCPKG_ROOT

1. 在 **"用户变量"** 区域，点击 **"新建"**
2. 变量名：`VCPKG_ROOT`
3. 变量值：`C:\vcpkg`
4. 点击 **"确定"**

### 步骤 7.3：添加 Ninja 到 PATH

1. 在 **"用户变量"** 中找到 `Path`，双击打开
2. 点击 **"新建"**
3. 输入：`C:\Tools\ninja`
4. 点击 **"确定"**

### 步骤 7.4：验证环境变量

1. **关闭所有命令提示符窗口**
2. 重新打开命令提示符
3. 验证：
   ```cmd
   echo %VCPKG_ROOT%
   ```
   应显示 `C:\vcpkg`

   ```cmd
   ninja --version
   ```
   应显示版本号

---

## 8. 克隆项目

### 步骤 8.1：选择项目位置

建议放在 `C:\Projects` 或 `D:\Projects`（路径不要有中文）

### 步骤 8.2：克隆代码

```cmd
cd C:\Projects
git clone <仓库地址> Flash
cd Flash
```

> 将 `<仓库地址>` 替换为实际的 Git 仓库 URL

---

## 9. 构建项目

### 步骤 9.1：打开开发者命令提示符

1. 按 `Win` 键，搜索 **"x64 Native Tools Command Prompt for VS"**
2. **右键** → **以管理员身份运行**

> **重要**：必须使用这个特殊的命令提示符，普通 cmd 无法编译！

### 步骤 9.2：进入项目目录

```cmd
cd C:\Projects\Flash
```

### 步骤 9.3：创建构建目录

```cmd
mkdir build
cd build
```

### 步骤 9.4：配置项目（CMake）

```cmd
cmake .. -G Ninja
```

**首次运行会自动安装依赖**，可能需要 5-15 分钟。

看到以下输出表示配置成功：
```
-- Configuring done
-- Generating done
-- Build files have been written to: C:/Projects/Flash/build
```

### 步骤 9.5：编译项目

```cmd
ninja
```

等待编译完成，看到类似输出即成功：
```
[xx/xx] Linking CXX executable RAND/stable_method_test.exe
```

---

## 10. 运行测试

在 build 目录下运行以下测试：

### 两相闪蒸测试

```cmd
.\RAND\2phase_tests.exe
```

预期输出：
```
[==========] Running 7 tests from 1 test suite.
...
[  PASSED  ] 7 tests.
```

### 三相闪蒸测试

```cmd
.\RAND\3phase_rand_test.exe
```

### 稳定法测试

```cmd
.\RAND\stable_method_test.exe
```

预期输出：
```
[  PASSED  ] 8 tests.
```

### 反应体系测试

```cmd
.\RAND\reactive_flash_test.exe
```

预期输出：
```
[  PASSED  ] 3 tests.
```

**所有测试通过即配置成功！**

---

## 11. 常见问题

### Q1: cmake 报错 "VCPKG_ROOT not set"

**原因**：环境变量未正确配置

**解决**：
1. 检查环境变量是否设置：`echo %VCPKG_ROOT%`
2. 如果为空，重新按第 7 节配置
3. **重新打开**命令提示符

---

### Q2: 找不到 "x64 Native Tools Command Prompt"

**原因**：Visual Studio Build Tools 未正确安装

**解决**：
1. 重新运行 `vs_BuildTools.exe`
2. 确保勾选 "使用 C++ 的桌面开发"
3. 重新安装

---

### Q3: ninja 命令找不到

**原因**：PATH 环境变量未配置

**解决**：
1. 确认 `C:\Tools\ninja\ninja.exe` 存在
2. 按第 7.3 节添加到 PATH
3. 重新打开命令提示符

---

### Q4: vcpkg 安装依赖失败

**解决**：手动安装依赖
```cmd
cd %VCPKG_ROOT%
.\vcpkg install eigen3:x64-windows
.\vcpkg install gtest:x64-windows
.\vcpkg install spdlog:x64-windows
.\vcpkg install openblas:x64-windows
```

---

### Q5: 运行测试时报错 "找不到 DLL"

**原因**：动态库未复制到可执行文件目录

**解决**：
构建系统会自动复制 DLL，如果仍有问题：
```cmd
copy ..\external\thermopack\lib\windows\*.dll .\RAND\
```

---

### Q6: 编译报错 "C++17 not supported"

**原因**：编译器版本过低

**解决**：
确保使用 Visual Studio 2019 或更高版本的 Build Tools

---

## 快速检查清单

配置完成后，确认以下命令都能正常执行：

```cmd
git --version          # Git 已安装
cmake --version        # CMake 已安装
ninja --version        # Ninja 已安装
echo %VCPKG_ROOT%      # 显示 C:\vcpkg
cl                     # 显示 MSVC 编译器信息
```

---

## 下一步

- 阅读 [CLAUDE.md](CLAUDE.md) 了解项目架构
- 查看 `RAND/tests/` 目录下的测试用例
- 参考 `RAND/include/rand_flash.hpp` 了解 API
