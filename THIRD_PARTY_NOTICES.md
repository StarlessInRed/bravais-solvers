# 来源与许可

PRISMS-PF core 由 PRISMS Center 上游维护，本仓不分发其完整源码；构建消费独立安装的 core。
应用接口及普通多晶粒 Allen–Cahn 方程依据 PRISMS-PF 应用方式实现。

`modules/grain_structure/upstream/` 的 `FloodFiller`、`OrderParameterRemapper`、
`SimplifiedGrainRepresentation` 头文件与实现来自 PRISMS-Center/phaseField 的
`150d85f74994eca7ea0d756d164e28eccef97b12`（v2.4.2），保持原算法及 LGPL-2.1-or-later 许可。
来源：https://github.com/PRISMS-Center/phaseField/tree/150d85f74994eca7ea0d756d164e28eccef97b12
`native.cpp` 和 `grain_growth_dream3d` 是面向 v4.0.1 生命周期及组织输入输出的我方适配；
构建核心仍为原版 v4.0.1，没有降级为 v2.4.2。

本仓应用和模块的我方修改由 wenjian 维护，Copyright (c) 2026 wenjian，
按各文件声明的 LGPL-2.1-or-later 提供。保留既有上游声明，许可全文见 LICENSE。
使用其他求解器时，其源码和许可由对应供给保留，本仓许可不覆盖它们。
