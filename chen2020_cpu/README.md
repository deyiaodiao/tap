# Chen 2020 CPU PBCD：RG_od 修复冻结基线

本目录复现 Chen、Liu、Zhang 与 Wang (2020) 的 CPU 并行路径型梯度投影
用户均衡配流算法（PBCD）：

> X. Chen, Z. Liu, K. Zhang, and Z. Wang,
> “A parallel computing approach to solve traffic assignment using path-based
> gradient projection algorithm,” *Transportation Research Part C*, 120,
> 102809, 2020.
> <https://doi.org/10.1016/j.trc.2020.102809>

当前版本是 2026-07-27 冻结的 PBCD 基线
`pbcd-rgod-baseline-v2`。它覆盖仓库中 2026-07-17 发布的旧基线。
旧版本仍可通过 Git 历史获得。

## 本次冻结修复了什么

旧基线使用“一次 GP 调流量”`maximumShift` 来决定 OD 是否继续进入
inner loop。论文 Fig. 2 写的是当前 OD 的相对间隙：

\[
RG_{od}=
\left|
1-\frac{q_{od}\,\pi_{od}}
{\sum_{k\in K_{od}} f_k c_k}
\right|.
\]

本版本先由当前路径费用与路径流计算 `RG_od`，仅当
`RG_od >= maxODGap` 时执行该 OD 的 GP 调流。默认筛选模式为
`relative-gap`。原来的 `flow-shift` 模式仅保留用于历史诊断，不属于冻结
基线协议。

算法的其余主要结构不变：

- 固定距离 OD 分块，默认每块 128 个 OD；
- 按起点并行生成最短路径；
- OD 私有流量增量与确定性合并；
- 并行刷新受影响链路费用和导数；
- 每次 inner loop 调流后更新链路流量与费用；
- 按起点并行计算全局 Relative Gap；
- Birmingham 只在程序内存中覆盖 BPR 参数为 `alpha=0.15, beta=4`，
  不修改 TNTP 原始文件。

## 冻结基线协议

冻结基线必须显式使用：

```text
--od-screening relative-gap
--gp-step 0.25
--od-per-block 128
--max-inner 1000
--full-check-frequency 100
--cost-scalar 60
```

源码为兼容既有调用仍保留历史 GP 默认值 0.30，但冻结基线不得依赖这个
隐式默认值。选择 0.25 的原因是它在 Winnipeg、Birmingham 和
Philadelphia 的 `1e-12` 回归中均收敛；0.30 在 Winnipeg 上曾达到
500 个 outer 上限但未收敛。

完整机器可读协议见
[`configs/pbcd_rgod_baseline_v2.json`](configs/pbcd_rgod_baseline_v2.json)。

## 环境与构建

需要 Linux、CMake 3.16 及以上、支持 C++17 和 OpenMP 的编译器。

```bash
cmake -S chen2020_cpu -B build/chen2020_cpu \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/chen2020_cpu -j6
ctest --test-dir build/chen2020_cpu --output-on-failure
```

程序位于 `build/chen2020_cpu/bin/chen2020_ue`。

`include/` 中的 6 个小写转发头文件用于兼容 Linux 区分大小写的文件
系统，不包含算法逻辑。

## 运行示例

Winnipeg：

```bash
./build/chen2020_cpu/bin/chen2020_ue pbcd \
  /path/to/wp Winnipeg /path/to/output/wp/solution \
  --threads 6 --convergence 1e-12 --max-outer 500 \
  --gp-step 0.25 --od-per-block 128 \
  --max-inner 1000 --full-check-frequency 100 \
  --cost-scalar 60 --od-screening relative-gap
```

Birmingham 还必须增加：

```text
--bpr-alpha 0.15 --bpr-beta 4
```

Philadelphia 不使用 BPR 覆盖。

## 冻结前高精度回归

以下结果均为 6 线程、目标 RG `1e-12`、GP 步长 0.25，并由独立验证器
从输出链路流重新计算：

| 网络 | Outer | Total inner | 独立 RG | Beckmann 目标函数 | 求解墙钟 |
|---|---:|---:|---:|---:|---:|
| Winnipeg | 41 | 7,434 | `7.9114462e-13` | `827911.4946299655` | 1.7435 s |
| Birmingham | 11 | 647 | `2.9718477e-13` | `193514.99208939876` | 22.5361 s |
| Philadelphia | 18 | 7,454 | `4.9988656e-13` | `271631357.31309846` | 225.0348 s |

三个网络均通过聚合链路流守恒验证。Birmingham 的两种步长解最大链路
流差为 `0.2268`。Philadelphia 的不同步长解在一条自由流时间为零的
连接边上可有较大流量差，但其费用恒为零，不影响 UE 目标与 RG。

详细验证范围、哈希和限制见
[`docs/VALIDATION.md`](docs/VALIDATION.md) 与
[`PROVENANCE.json`](PROVENANCE.json)。

## 验证边界

- 独立验证覆盖链路流守恒、最短路、Relative Gap 和 Beckmann 目标。
- Winnipeg 与 Philadelphia 本轮未保留完整 OD 路径流文件，因此不能
  表述为“每个 OD 的路径流已独立审计”。
- Greedy/TNM 是研究中指定的独立串行控制，不是论文作者公开的 iGP
  源码。
- 论文作者未公开原始实现，因此本仓库验证的是论文公式、收敛行为和
  UE 条件，不能声明与作者代码逐行一致。
- 网络数据和大规模原始结果不随本仓库分发。

## 权利说明

本目录包含不同来源的研究代码，当前未附统一开源许可证。详情见
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)。
