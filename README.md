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

## 项目结构

```
Flash/
├── Thermo/              # 热力学后端层
├── PhaseStability/      # 相稳定性分析
├── RAND/                # 多相闪蒸求解器
└── external/thermopack/ # 内嵌ThermoPack库
```

## 架构设计

```
RAND (闪蒸求解器)
    ↓
PhaseStability (稳定性分析)
    ↓
Thermo (热力学后端)
    ↓
外部库 (ThermoPack, Eigen3)
```

## 文档

- [使用指南](GUIDE.md) - 环境配置和快速开始
- [开发文档](CLAUDE.md) - 详细的开发者文档

## 许可证

[在此添加许可证信息]
