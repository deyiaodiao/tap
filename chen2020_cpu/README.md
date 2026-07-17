# Chen 2020 CPU PBCD 复现

本目录复现 Chen、Liu、Zhang 与 Wang (2020) 提出的 CPU 并行路径型梯度投影用户均衡配流方法：Parallel Block-Coordinate Descent（PBCD）。论文为 [*A parallel computing approach to solve traffic assignment using path-based gradient projection algorithm*](https://doi.org/10.1016/j.trc.2020.102809)。

当前代码是 2026-07-17 冻结的算法基线。算法源文件与冻结版本逐字节一致；发布副本未包含本机构建产物、网络数据和大规模原始输出。身份信息见 [`PROVENANCE.json`](PROVENANCE.json)。

`include/` 中另外加入了 6 个小写转发头文件，使旧 TNM/NEWRAN 大写文件名在区分大小写的 Linux 文件系统上也能构建；这些文件不包含算法逻辑，且没有改动冻结头文件。

## 主要实现

- 固定距离 OD 分块，默认每块 128 个 OD；
- 按起点并行生成最短路径；
- OD 私有流量增量与确定性合并；
- 并行刷新受影响链路费用和导数；
- 按起点并行计算 Relative Gap；
- Birmingham 网络加载后可在内存中统一覆盖 BPR 参数为 `alpha=0.15, beta=4`，不修改原始数据。

## 环境与构建

需要 Linux、CMake 3.16 及以上、支持 C++17 和 OpenMP 的编译器。已验证环境为 WSL2 Ubuntu 24.04、GCC 13.3、6 个物理核。

从仓库根目录执行：

```bash
cmake -S chen2020_cpu -B build/chen2020_cpu -DCMAKE_BUILD_TYPE=Release
cmake --build build/chen2020_cpu -j6
ctest --test-dir build/chen2020_cpu --output-on-failure
```

程序位于 `build/chen2020_cpu/bin/chen2020_ue`。

## 运行示例

输入目录中应包含 `<network-name>_net.tntp` 和 `<network-name>_trips.tntp`。输出前缀必须位于输入目录之外。

Birmingham：

```bash
./build/chen2020_cpu/bin/chen2020_ue pbcd \
  /path/to/Birmingham-England Birmingham /path/to/output/bh/solution \
  --threads 6 --convergence 1e-12 \
  --gp-step 0.3 --od-per-block 128 \
  --max-inner 1000 --full-check-frequency 100 \
  --cost-scalar 60 \
  --bpr-alpha 0.15 --bpr-beta 4
```

Philadelphia 不使用统一 BPR 覆盖，因此省略最后两个参数。性能实验把 `--convergence` 改为 `1e-5`；高精度实验使用 `1e-12`。完整冻结协议见 [`configs/chen2020_formal_v1.json`](configs/chen2020_formal_v1.json)。

## 已验证结果

正式 v1 协议在 6 个物理核上采用一次预热和五次计时，时间口径为报告写出前的墙钟中位数：

| 网络 | PBCD 1 线程 | PBCD 6 线程 | PBCD 加速比 | Greedy 串行基准 | Greedy/PBCD6 |
|---|---:|---:|---:|---:|---:|
| Birmingham | 71.600 s | 14.832 s | 4.8276x | 69.173 s | 4.6639x |
| Philadelphia | 2148.245 s | 533.985 s | 4.0230x | 1102.053 s | 2.0638x |

高精度独立链路验证结果：

| 网络 | 6 线程 Relative Gap | Beckmann 目标值 | 验证范围 |
|---|---:|---:|---|
| Birmingham | `6.9083e-15` | `193514.9920894` | 链路守恒；另完成全部 470,805 个 OD 的路径流审计 |
| Philadelphia | `8.4477e-14` | `271631357.313099` | 链路守恒与独立最短路/RG 重算 |

结果解释、限制和验证边界见 [`docs/VALIDATION.md`](docs/VALIDATION.md)。

## 重要边界

- `PBCD 1 线程` 是同一分块/Jacobi 实现的单线程控制，不等同于论文中的逐 OD iGP；正式并行加速比以它为分母。
- Greedy/TNM 是本次复现指定的独立串行基准，不应表述为论文原始 iGP 源码。
- 论文作者未公开实现，因此这里只能验证公式、收敛轨迹、UE 条件和性能区间，不能声明代码级完全一致。
- 网络数据不在本仓库中；Transportation Networks 数据集要求用于学术研究并标注来源。

## 引用

```text
X. Chen, Z. Liu, K. Zhang, and Z. Wang,
"A parallel computing approach to solve traffic assignment using
path-based gradient projection algorithm," Transportation Research
Part C: Emerging Technologies, vol. 120, 102809, 2020.
https://doi.org/10.1016/j.trc.2020.102809
```
