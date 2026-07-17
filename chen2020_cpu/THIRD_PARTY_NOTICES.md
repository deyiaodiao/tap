# 第三方来源与权利说明

本目录不是从零编写的独立代码库。PBCD 实现建立在既有 TNM/Greedy 交通配流框架的必要源文件之上。发布者不对第三方文件主张所有权。

## TNM（Transportation Network Modelling）框架

`include/TNM_Header.h` 的原始文件头标注：

- 项目：TNM — Transportation network modelling
- 作者：Yu Nie
- 机构：UC Davis
- 最后更新：2004 年 9 月

当前取得的 TNM 源文件中没有发现明确的统一许可证文本。`include/My_Predicate.h` 还注明其基于一个来自网页、但原作者姓名未被记录的程序。公开可见不等于获得开源许可；在复制、再分发或商业使用前，应由使用者进一步确认相应权利。

## Newran

`include/randlib/` 包含 Robert B. Davies 的 NEWRAN02B 头文件（2002-07-22）。作者的 [Newran 文档](https://www.robertnz.net/nr02doc.htm) 声明该库使用不受限制、不承担责任，并要求分发源代码时明确哪些部分属于作者、可从互联网免费取得。

本仓库据此保留 Newran 文件的原始标识，并在此明确归属。Newran 文档同时说明其中部分数学和排序代码改编自其他出版物；如扩展或重新分发完整 Newran 实现，应继续核查这些来源。

## 本复现新增部分

`src/PBCD_Algorithm.cpp`、`include/PBCD_Algorithm.h`、`src/driver.cpp`、对应测试及发布文档属于本次 Chen 2020 复现工作形成的新增或适配内容。由于它们与未明确许可的 TNM 框架共同构建，仓库当前不附统一开源许可证。

此文件是来源记录，不构成法律意见或对任何权利的保证。
