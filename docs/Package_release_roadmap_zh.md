# BoundaryTMatrixMethods.jl 正式成包路线

## 1. 包的定位

建议把包定位成：

> 一个把 SVM/PMM 边界匹配算法转化为 `TransitionMatrices.jl` T 矩阵对象的 Julia 扩展包，并提供 residual diagnostics、PMM gate 和 IITM fallback。

不要把它写成“发明了 SVM/PMM”。SVM 和 PMM 本身已有文献。你的贡献应该写成：

- Julia 实现；
- 对接 `TransitionMatrices.jl`；
- 输出标准 T 矩阵对象；
- 可以直接算振幅矩阵、Mueller 矩阵、散射/消光截面；
- 提供 PMM 是否可信的自动诊断；
- 给 IITM 提供可控的加速/辅助入口。

## 2. 包名建议

当前包名：

```text
BoundaryTMatrixMethods.jl
```

这个名字是比较稳的，含义清楚，不会过度声明。

备选名字：

```text
BoundaryMatchingTMatrix.jl
SVMPMMTMatrix.jl
PointMatchingTMatrix.jl
BoundaryTMatrix.jl
```

我更建议保留 `BoundaryTMatrixMethods.jl`。原因：

- 不只包含 PMM，也包含 SVM；
- 不只求解，也包含 residual diagnostics 和 selector；
- 后面还能扩展到更多 boundary-matching 方法。

## 3. 包结构

当前结构已经接近正式包：

```text
BoundaryTMatrixMethods.jl/
  Project.toml
  README.md
  src/
    BoundaryTMatrixMethods.jl
    boundary_system.jl
    solvers.jl
    diagnostics.jl
    correction.jl
    selector.jl
    metrics.jl
  test/
    runtests.jl
    iitm_integration.jl
  examples/
    basic_validation.jl
    pmm_gate_demo.jl
  docs/
  scripts/
  paper_validation/
```

正式发布前建议拆成两层：

### 包核心仓库

保留：

- `src/`
- `test/`
- `examples/`
- `docs/`
- `README.md`
- `Project.toml`

### 论文验证材料

可以保留在同一个仓库，也可以单独放：

- `scripts/`
- `paper_validation/data/`
- `paper_validation/figures/`
- `paper_validation/reports/`

如果投 CPC，最好保留验证脚本和机器可读数据，因为这体现 reproducibility。

## 4. API 设计

正式包最关键的是 API 稳定。

建议主 API 保持简单：

```julia
Tsvm = svm_tmatrix(shape, lambda, nmax, ngauss)
Tpmm = pmm_tmatrix(shape, lambda, nmax, ngauss)
```

这两个函数必须返回：

```julia
TransitionMatrices.AxisymmetricTransitionMatrix
```

这样用户可以直接：

```julia
S = amplitude_matrix(Tpmm, theta_i, phi_i, theta_s, phi_s; λ=lambda)
F = scattering_matrix(Tpmm, lambda, angles)
Csca = scattering_cross_section(Tpmm, lambda)
Cext = extinction_cross_section(Tpmm, lambda)
```

诊断 API：

```julia
gate = diagnose_pmm_applicability(shape, lambda, nmax; ngauss=120)
```

自动选择 API：

```julia
result = auto_tmatrix(shape, lambda, nmax; Nr=6, Ntheta=40)
T = result.T
```

这里 `result.method` 应该告诉用户最终用了 `:pmm` 还是 `:iitm`。

## 5. 对接 TransitionMatrices.jl 的原则

你的包不要重新定义 T 矩阵类型。

应该依赖：

```julia
TransitionMatrices.AxisymmetricTransitionMatrix
```

也就是说，SVM/PMM 只负责生成 block：

```text
m = 0, 1, 2, ..., nmax
```

然后组装成 `AxisymmetricTransitionMatrix`。

这样你天然继承：

- `amplitude_matrix`
- `scattering_matrix`
- `phase_matrix`
- `scattering_cross_section`
- `extinction_cross_section`
- `orientation_average`
- `transition_matrix_iitm`

这就是包的最大价值。

## 6. 必须写清楚适用范围

README 和论文里都要明确：

当前版本适用：

- axisymmetric particles；
- homogeneous particles；
- shapes 支持 `TransitionMatrices.gaussquad(shape, ngauss)`；
- shape 需要有 `shape.m` 作为相对复折射率；
- spheroid、Chebyshev 已验证。

当前版本不应过度承诺：

- 非轴对称粒子；
- 多层粒子；
- 大尺寸参数强形变粒子；
- PMM 在所有情况下都稳定；
- boundary residual 等于最终光学量误差。

## 7. 必须有 gate

这不是附加功能，而是论文/包的核心卖点。

PMM/SVM 在小尺寸或温和形变下可以很好，但大尺寸参数时会失效。例如：

- `a=0.5, c=0.8, λ=1`：PMM/SVM 截面还不错，但振幅矩阵误差约 6.9%；
- `a=0.5, c=0.8, λ=0.1`：PMM/SVM 高阶直接失稳，PMM 出现 `Csca >> Cext`。

所以包里必须保留：

```julia
diagnose_pmm_applicability
auto_tmatrix
```

这能把你的包从“算法翻译”变成“可靠性工作流”。

## 8. 测试体系

最低测试集：

1. 球体
   - SVM/PMM 与 EBCM/Mie 达到机器精度或接近机器精度。

2. 温和椭球
   - Mueller 矩阵误差小；
   - Csca/Cext 合理；
   - gate 接受。

3. Chebyshev n=2
   - gate 接受；
   - Mueller 曲线与 EBCM 接近。

4. stress cases
   - 强椭球；
   - Chebyshev n=4；
   - gate 应拒绝 PMM；
   - `auto_tmatrix` 应回退 IITM。

5. 振幅矩阵测试
   - 至少一个 `amplitude_matrix` 对照 EBCM；
   - 不只看 `scattering_matrix`。

## 9. 文档结构

README 建议包含：

- installation；
- quick start；
- SVM/PMM examples；
- PMM gate example；
- 与 `TransitionMatrices.jl` 对接说明；
- failure mode；
- citation。

docs 建议包含：

- `PMM_residual_gate.md`
- `CPC_article_storyline_zh.md`
- `CPC_manuscript_blueprint.md`
- `Package_release_roadmap_zh.md`

后面正式注册或投稿时，可以再加：

- `docs/src/index.md`
- `docs/make.jl`
- Documenter.jl 文档站。

## 10. CPC 文章怎么讲

文章标题可以是：

> BoundaryTMatrixMethods.jl: SVM/PMM boundary-matching T-matrix generation with residual diagnostics and IITM fallback in Julia

核心贡献：

```text
SVM/PMM algorithms -> TransitionMatrices.jl T matrix -> validation -> residual gate -> IITM fallback
```

文章不要主打：

```text
我们发明了 PMM/SVM
```

要主打：

```text
我们把 PMM/SVM 做成了可复现、可诊断、可接入 IITM 的 Julia T 矩阵包
```

## 11. 正式发布前的 TODO

高优先级：

- 把 `TransitionMatrices.jl` 依赖从本地路径改成正式可解析依赖；
- 补充 `compat`，例如 `TransitionMatrices = "...版本..."`；
- 确认 `Pkg.test()` 在干净环境通过；
- 添加 CI；
- README 加完整安装说明；
- 给每个 public function 加 docstring；
- 增加 `amplitude_matrix` 验证案例；
- 做参数扫图：aspect ratio、wavelength、size parameter、Chebyshev amplitude。

中优先级：

- 增加 Documenter.jl 文档；
- 增加 benchmark；
- 增加 Zenodo DOI；
- 增加 CITATION.cff；
- 增加 LICENSE；
- 增加 CPC `PROGRAM SUMMARY`。

## 12. 最终建议

你的做法应该是：

1. 继续用 `BoundaryTMatrixMethods.jl` 作为包名；
2. 核心算法只负责输出 `AxisymmetricTransitionMatrix`；
3. 所有后处理都调用 `TransitionMatrices.jl`；
4. 把 residual gate 作为包的安全机制；
5. 论文主线写“软件接口 + 可靠性诊断 + 自动 fallback”；
6. 不要承诺 PMM/SVM 对所有大尺寸粒子都有效。

这样这个包才有机会从“代码复现”变成“可发表的软件贡献”。
