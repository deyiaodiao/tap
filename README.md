# TAP

交通配流（Traffic Assignment Problem）算法复现与后续研究成果仓库。

## 项目

- [`chen2020_cpu/`](chen2020_cpu/)：Chen et al. (2020) 路径型梯度投影
  用户均衡配流算法的 CPU 并行 PBCD 复现。当前冻结基线为
  `pbcd-rgod-baseline-v2`（2026-07-27），OD 筛选采用论文所写的
  `RG_od >= maxODGap`。

旧 PBCD 基线已经由当前版本覆盖，但仍可通过 Git 历史追溯。后续成果
将作为与 `chen2020_cpu/` 并列的一级目录加入本仓库。

## 数据

仓库不重复分发大型路网和实验原始输出。Birmingham、Chicago Sketch、Philadelphia 与 Winnipeg 的 TNTP 数据可从 [Transportation Networks for Research](https://github.com/bstabler/TransportationNetworks) 获取；使用数据时请遵守其学术研究和来源标注要求。

## 权利说明

本仓库包含具有不同来源的研究代码，当前未附统一开源许可证。具体来源和已知权利声明见 [`chen2020_cpu/THIRD_PARTY_NOTICES.md`](chen2020_cpu/THIRD_PARTY_NOTICES.md)。
