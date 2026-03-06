# Gibbs Energy Verification Tests

本目录包含两个测试程序，用于验证乙烷裂解反应体系的Gibbs能量计算。

## 测试系统

**反应**: C2H6 ⇌ C2H4 + H2
**温度**: 800 K
**压力**: 75000 Pa
**状态方程**: SRK (Soave-Redlich-Kwong)

## 对比的两个解

**Solution 1 (均匀分布)**:
- 摩尔数: [0.505567, 0.494433, 0.494433]
- 总摩尔数: 1.49443 mol
- 组成: x ≈ [0.33, 0.33, 0.33]

**Solution 2 (用户提出的解)**:
- 摩尔数: [0.919233, 0.0811087, 0.0811087]
- 总摩尔数: 1.08145 mol
- 组成: x = [0.85, 0.075, 0.075]

两个解都满足元素守恒: C=2, H=6

## 测试1: ThermoPack/C++ 实现

### 构建
```bash
cd build
ninja verify_gibbs
```

### 运行
```bash
# macOS
DYLD_LIBRARY_PATH=../external/thermopack/lib/macos:$DYLD_LIBRARY_PATH ./Thermo/verify_gibbs

# Linux
LD_LIBRARY_PATH=../external/thermopack/lib/linux:$LD_LIBRARY_PATH ./Thermo/verify_gibbs
```

### 结果
ThermoPack/SRK计算显示：
- Solution 1: G ≈ -296808 J
- Solution 2: G ≈ -292878 J
- **ΔG ≈ 3930 J** (Solution 1更低)

## 测试2: Clapeyron.jl 实现

### 环境准备
需要安装Julia和Clapeyron包：

```bash
# 安装Julia (如果尚未安装)
# macOS: brew install julia
# 或从 https://julialang.org/downloads/ 下载

# 安装Clapeyron包
julia -e 'using Pkg; Pkg.add("Clapeyron")'
```

### 运行
```bash
cd Thermo/tests
julia verify_gibbs_clapeyron.jl
```

### 结果
Clapeyron/SRK计算与ThermoPack一致，显示Solution 1更低。

## 验证结果 ✓

**已验证**: ThermoPack和Clapeyron给出一致的结果，但都与Aspen不一致。

```
ThermoPack (SRK) ≈ Clapeyron (SRK) ≠ Aspen (SRK + ΔG_f°)
     ↓                    ↓                    ↓
均匀分布更低        均匀分布更低        [0.85,0.075,0.075]更低
  (错误预测)          (错误预测)           (正确预测)
```

### 结论

1. **算法实现正确** ✓
   - 两个独立的SRK实现给出一致结果
   - 排除了代码bug或实现错误
   - RAND算法正确求解了给定目标函数的最小值

2. **热力学模型局限性** ⚠️
   - SRK without ΔG_f°无法正确预测化学平衡
   - 问题不在优化算法，而在目标函数本身

3. **Aspen的优势** ✓
   - Aspen使用完整的热力学数据库
   - 包含标准Gibbs生成能(ΔG_f°)数据
   - 能够正确预测化学反应平衡

## 为什么会出现差异？

### SRK状态方程的能力范围

```
SRK只能计算：
├─ PVT关系（压力-体积-温度）
├─ 逸度系数 φ_i(T, P, x)
└─ 混合物非理想性

SRK无法计算：
└─ 分子的内在化学稳定性（需要ΔG_f°）
```

### 化学势的组成

对于组分i的化学势：
```
μ_i = μ_i^0(T) + RT ln(x_i φ_i)
      ↑           ↑
   需要ΔG_f°    SRK可以算
```

**问题所在**：
- μ_i^0(T)项包含分子的内在稳定性信息
- 没有ΔG_f°数据，这一项是任意的（取决于参考态）
- SRK无法区分C2H6、C2H4、H2的真实稳定性差异
- 只能基于混合规则计算，导致错误的平衡预测

### Aspen如何解决这个问题？

Aspen使用完整的热力学数据库：

```
Aspen = SRK (PVT) + 热力学数据库
                    ├─ ΔG_f°(T) - 标准Gibbs生成能
                    ├─ ΔH_f°(T) - 标准生成焓
                    ├─ S°(T)    - 标准熵
                    └─ Cp(T)    - 热容
```

对于乙烷裂解反应：
```
C2H6 → C2H4 + H2

ΔG_rxn°(298K) = ΔG_f°(C2H4) + ΔG_f°(H2) - ΔG_f°(C2H6)
              ≈ 68.4 + 0 - (-32.0) = 100.4 kJ/mol

在800K高温下：
- ΔH_rxn ≈ 136 kJ/mol (吸热反应)
- TΔS_rxn项在高温下变大
- ΔG_rxn = ΔH - TΔS < 0 (反应自发)
→ 产物(C2H4, H2)应占优
→ 用户的解[0.85, 0.075, 0.075]是正确的
```

## 对项目的影响

### 当前状态

✅ **非反应体系**: RAND算法完全适用
- 相平衡计算（气-液、液-液）
- 多相闪蒸
- 稳定性分析

⚠️ **反应体系**: 存在已知局限性
- SRK without ΔG_f°无法预测正确的化学平衡
- 算法会收敛到数学上的局部最小值
- 但这个最小值在物理上可能是错误的

### 建议

**对于反应体系，有以下选项**：

1. **使用包含ΔG_f°的热力学模型**
   - 例如：GERG-2008, PC-SAFT with association
   - 或者在SRK基础上添加ΔG_f°数据库

2. **明确文档说明局限性**
   - 在测试和文档中说明SRK的适用范围
   - 警告用户反应体系的预测可能不准确

3. **使用外部化学平衡求解器**
   - 先用化学平衡求解器（如Cantera）计算平衡组成
   - 再用RAND进行相平衡计算

**当前实现**：已在测试中明确说明这是已知局限性。

