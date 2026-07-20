# BoundaryTMatrixMethods.jl 文章内容规划（中文）

## 文章定位

这篇文章不要写成“SVM 和 PMM 新算法”。这两个方法本身不是新的，审稿人很容易指出已有文献。更稳的定位是：

> 我们做了一个可复现的 Julia 软件包，把 SVM/PMM 边界匹配算法变成 `TransitionMatrices.jl` 可直接使用的 T 矩阵生成器，并且加入 residual diagnostic gate，自动判断 PMM 什么时候可信、什么时候应该回退到 IITM。

也就是说，文章的真正主线是：

```text
算法实现 -> T 矩阵接口 -> Mueller/cross-section 验证 -> PMM 可靠性诊断 -> 自动 selector -> 可复现软件包
```

这个叙事比“我把 scattpy 改成 Julia”强很多，也比“PMM 一定加速 IITM”稳很多。

## 建议题目

首选：

> BoundaryTMatrixMethods.jl: boundary-matching T-matrix generation with residual diagnostics and IITM fallback in Julia

备选：

> A Julia implementation of SVM/PMM T-matrix methods with no-reference reliability diagnostics for axisymmetric light scattering

更偏 CPC 软件风格：

> BoundaryTMatrixMethods.jl: a reproducible SVM/PMM T-matrix interface for TransitionMatrices.jl

## 摘要要讲什么

摘要应该包含五句话：

1. T 矩阵方法在非球形粒子散射中很重要，但不同算法之间的软件接口不统一。
2. 本文提出 `BoundaryTMatrixMethods.jl`，实现 SVM/PMM 边界匹配 T 矩阵生成。
3. 该包直接返回 `TransitionMatrices.jl` 的 `AxisymmetricTransitionMatrix`，可以无缝计算 Mueller 矩阵和截面。
4. 为避免 PMM 被盲目使用，本文提出 no-reference residual gate，结合 refined-grid stability、boundary residual、projected residual 和 passivity。
5. 球体、温和椭球、Chebyshev 粒子的验证显示，gate 接受的 PMM/SVM 结果可靠；强形变 stress cases 会被 gate 拒绝并自动回退到 IITM。

## 文章贡献点

建议明确列 5 个贡献：

1. Julia 中实现 SVM 和 PMM 边界匹配 T 矩阵算法。
2. 与 `TransitionMatrices.jl` 的 T 矩阵类型衔接。
3. 提供六个独立 Mueller 矩阵元素、散射截面、消光截面的完整验证流程。
4. 提出 PMM residual diagnostic gate，说明 PMM 什么时候可以用、什么时候不能用。
5. 提供自动 selector：gate 通过时用 PMM，否则回退到 IITM。

## 正文章节安排

### 1. Introduction

这一节回答“为什么值得做”。

要讲：

- 非球形粒子散射需要 T 矩阵。
- EBCM、IITM、SVM、PMM 都能在不同场景中生成 T 矩阵，但软件接口常常分散。
- 用户真正需要的是一个统一 T 矩阵对象，然后可以继续算 Mueller matrix、cross sections、orientation averaging 等。
- PMM/SVM 的问题不是“有没有算法”，而是“能否可靠接入现有 T 矩阵生态，并知道什么时候结果可信”。

建议最后一句：

> The main contribution of this work is therefore not a new boundary-matching algorithm, but a reproducible software layer that couples SVM/PMM T-matrix generation to residual diagnostics and IITM fallback.

对应图：

- `nature_fig01_graphical_abstract`
- `nature_fig10_manuscript_roadmap`

### 2. Mathematical Formulation

这一节回答“算法如何形成 T 矩阵”。

要讲：

- 入射场、散射场、内部场都用 vector spherical wave functions 展开。
- T 矩阵把入射展开系数映射到散射展开系数。
- 轴对称粒子可以按 azimuthal order m 分块处理。
- 边界条件是切向电场和磁场连续。
- SVM 可以理解为在边界上建立 surface variable system。
- PMM 可以理解为用边界采样点/加权残差来最小化边界不连续。
- 最终每个 m block 被组装成 `AxisymmetricTransitionMatrix`。

这里不需要推太多公式，重点是和软件接口对齐。

对应图：

- `nature_fig01_graphical_abstract`
- `nature_fig09_api_reproducibility`

### 3. Software Design

这一节回答“包怎么用”。

要讲：

- `svm_tmatrix(shape, lambda, nmax, ngauss)`
- `pmm_tmatrix(shape, lambda, nmax, ngauss)`
- `diagnose_pmm_applicability(...)`
- `auto_tmatrix(...)`
- `boundary_residual(...)`
- `projected_boundary_residual(...)`
- `mueller_relative_error(...)`

要强调：

- 输出是 `TransitionMatrices.jl` 原生 T 矩阵对象。
- 因此可以直接调用 `scattering_matrix`、`scattering_cross_section`、`extinction_cross_section`。
- 数据和图都是脚本可复现的。

对应图：

- `nature_fig09_api_reproducibility`

### 4. Residual Diagnostic Gate

这是文章最关键的创新点之一。

要讲清楚：

PMM 有时候很好，有时候会错。不能只看 boundary residual，因为低 residual 只能说明当前离散系统解得一致，不保证 T 矩阵已经物理收敛。

所以 gate 使用四个指标：

1. `convergence_error`
   - 基础 PMM 与 refined nmax/ngauss PMM 的 T 矩阵差异。
   - 这是最重要的 no-reference 指标。

2. `boundary_residual`
   - 原始边界条件残差。

3. `projected_boundary_residual`
   - 投影掉内部场自由度后的残差。

4. `energy_violation`
   - 被动粒子一致性检查，即是否出现明显 `Csca > Cext`。

gate 逻辑：

```text
if all diagnostics pass:
    use PMM
else:
    fall back to IITM
```

对应图：

- `nature_fig06_gate_diagnostics`

### 5. Validation Setup

这一节回答“拿什么证明它对”。

当前验证集：

- sphere
- mild prolate spheroid
- mild oblate spheroid
- Chebyshev n=2
- strong prolate spheroid
- strong oblate spheroid
- Chebyshev n=4 positive
- Chebyshev n=4 negative

参考方法：

- EBCM 作为主要 reference。
- IITM coarse/fine 用作集成测试和 fallback reference。

验证指标：

- T-matrix relative error
- Mueller-matrix relative error
- F11 relative error
- Csca relative error
- Cext relative error
- boundary residual
- projected residual
- gate status

对应图：

- `nature_fig02_case_library`
- `nature_fig03_accuracy_forest`

### 6. Results

这一节建议分成四个结果小节。

#### 6.1 Sphere sanity check

讲：

- 球体是最基本验证。
- SVM/PMM 与参考结果达到机器精度量级。
- 说明 T 矩阵 block 组装、归一化、Mueller 接口没有明显错误。

可引用数值：

- sphere PMM Mueller error: `7.350e-16`
- sphere SVM Mueller error: `1.020e-15`

对应图：

- `nature_fig03_accuracy_forest`
- supplementary sphere six-element figure

#### 6.2 Gate-accepted non-spherical cases

讲：

- mild prolate、mild oblate、Chebyshev n=2 被 gate 接受。
- 六个独立 Mueller 元素和 EBCM 基本重合。
- 这证明 SVM/PMM 生成的 T 矩阵在这些情形下可以直接用于物理观测量。

可引用数值：

- mild prolate PMM Mueller error: `7.935e-04`
- mild oblate PMM Mueller error: `1.405e-03`
- Chebyshev n=2 PMM Mueller error: `8.680e-05`

对应图：

- `nature_fig04_mueller_accepted_atlas`
- `nature_fig08_observable_decomposition`

#### 6.3 Stress cases and PMM failure

讲：

- strong prolate、strong oblate、Chebyshev n=4 是 stress cases。
- PMM 的 boundary residual 可能不太大，但 refined-grid stability 失败。
- 因此 PMM 不能盲目使用。
- selector 在这些 cases 自动回退 IITM。

可引用数值：

- strong prolate PMM Mueller error: `1.184e-02`
- strong oblate PMM Mueller error: `1.618e-02`
- Chebyshev n=4 positive PMM Mueller error: `8.299e-04`
- Chebyshev n=4 positive IITM coarse error: `6.129e-05`

对应图：

- `nature_fig05_stress_deviation_atlas`
- `nature_fig06_gate_diagnostics`

#### 6.4 Convergence and observables

讲：

- 随 nmax 增加，accepted cases 的 SVM/PMM Mueller error 呈下降或稳定趋势。
- T 矩阵元素误差、Mueller 矩阵误差和截面误差并不完全等价。
- 因此文章报告多种 observables 是必要的。

对应图：

- `nature_fig07_convergence`
- `nature_fig08_observable_decomposition`

### 7. Discussion

这一节要主动讲限制，避免审稿人抓住说过度宣传。

限制：

- 当前实现主要面向 axisymmetric homogeneous particles。
- 强形变或高阶 Chebyshev 情况下 PMM 可能不稳定。
- gate 是 no-reference diagnostic，不是数学严格误差界。
- 当前 runtime 图只能作为 snapshot，正式投稿前需要 warm-start benchmark。
- 还需要更大参数扫图：size parameter、aspect ratio、Chebyshev amplitude、refractive index。

但要强调：

- 正因为 PMM 会失败，所以 residual gate 和 IITM fallback 是必要贡献。
- 包的价值不是“PMM 永远更准”，而是“PMM 可信时自动用，不可信时自动停”。

### 8. Conclusion

结论可以这样收：

> BoundaryTMatrixMethods.jl provides a practical bridge between SVM/PMM boundary matching and the Julia T-matrix ecosystem. By combining TransitionMatrices-compatible outputs with a no-reference residual gate and automatic IITM fallback, the package makes PMM/SVM T matrices usable in a controlled and reproducible way.

## 主图安排

建议主文 6-8 张图，不一定全部塞进主文。Nature-style 图组里 10 张图的用途如下：

1. `nature_fig01_graphical_abstract`
   - 图形摘要/Introduction。

2. `nature_fig02_case_library`
   - Validation setup。

3. `nature_fig03_accuracy_forest`
   - 总体精度图。

4. `nature_fig04_mueller_accepted_atlas`
   - 主正确性图，必须放主文。

5. `nature_fig05_stress_deviation_atlas`
   - PMM failure mode，建议主文。

6. `nature_fig06_gate_diagnostics`
   - residual gate 的核心图，必须放主文。

7. `nature_fig07_convergence`
   - 可放主文或补充。

8. `nature_fig08_observable_decomposition`
   - 可放主文，证明验证不只看 Mueller。

9. `nature_fig09_api_reproducibility`
   - 软件设计图，可放 Software Design。

10. `nature_fig10_manuscript_roadmap`
    - 内部规划图，不一定放论文。

## 表格安排

建议至少有 4 张表：

1. Algorithm/API table
   - function
   - role
   - output
   - uses EBCM reference or not

2. Benchmark particle table
   - shape
   - parameters
   - nmax
   - quadrature
   - group: accepted/stress

3. Accuracy summary table
   - SVM Mueller error
   - PMM Mueller error
   - T error
   - Csca/Cext error

4. Gate diagnostics table
   - convergence error
   - boundary residual
   - projected residual
   - energy violation
   - selected method

## 现在还需要补什么才能更像 CPC 投稿

最高优先级：

1. 参数扫图
   - aspect ratio sweep
   - Chebyshev amplitude sweep
   - refractive index sweep
   - size parameter sweep

2. runtime benchmark
   - 去掉 Julia 编译影响
   - warm-up 后重复多次
   - 报 median 和 interquartile range

3. README 和 examples 完善
   - 一键运行案例
   - 一键重画图
   - 测试命令

4. 与已有开源实现比较
   - 至少说明已有开源软件做了什么，本包补了什么。

5. residual gate 阈值校准
   - 当前阈值是 conservative validation。
   - 更强的文章需要二维或三维 parameter map 来说明阈值如何区分 safe/unsafe。

## 投稿机会判断

如果只写“SVM/PMM Julia 版本”，机会不大。

如果写成：

```text
Julia T-matrix generator + TransitionMatrices interface + residual diagnostic gate + auto IITM fallback + reproducible validation package
```

就有 CPC 软件论文的机会。

当前材料大概是可投稿雏形，但还不够“稳投”。补完参数扫图和正式 benchmark 后，机会会明显提高。
