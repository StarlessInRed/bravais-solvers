# Bravais Solvers

Bravais 的可独立编译数值源码。当前提供 PRISMS-PF 普通曲率驱动晶粒长大应用
`grain_growth_dream3d`，支持二维/三维、一阶有限元、运行时 count 以及上游晶粒分配和通道重映射。
晶粒编号、晶体取向与序参量通道不是同一个量；多个充分分离的晶粒可以复用通道。

可选 `grain_growth_gnd` 提供二维 GND proxy 驱动再结晶：母晶逐网格储能、初始晶核和逐帧晶粒身份。
它使用同一 PRISMS-PF core，编译时将上面应用路径替换为 `applications/prismspf/grain_growth_gnd`，
运行产物为 `bin/grain_growth_gnd`。产品负责生成 `initial_fields.txt`、`initial_identity.txt`、
`stored_energy.txt` 和 PRM；本切片固定无回复、无演化中形核，RX 密度为零，不包含完整缺陷模型。
GUI 中的“GND再结晶”使用二维 EBSD 的 GND proxy，不表示完整三维位错测量。

## 构建与输入

需要独立安装的 PRISMS-PF 4.0.1（提交 `62b31b703734e43235d595525e9aecd58b6bdf64`）、
兼容的 deal.II、C++20 编译器和 CMake 3.25 以上。

```sh
cmake -S applications/prismspf/grain_growth_dream3d -B build -DPRISMS-PF_DIR=<core>/lib/cmake/PRISMS-PF
cmake --build build --parallel 2
cmake --install build --prefix <install>
```

在含 `parameters.prm` 的工作目录运行 `<install>/bin/grain_growth_dream3d`。
输入由调用产品生成：PRM、组织晶粒编号 `feature_ids.txt`；原生程序负责通道分配和平滑。
count 由输入决定，不需要修改 C++。输出包括 VTU 连续场和逐帧晶粒身份。
本应用不接受储能、回复或出生输入；运行成功不代替材料标定或网格独立性检验。

`modules/grain_growth/model.*` 实现普通长大势阱；`modules/grain_structure` 接入上游分配、换场和 IO。
输入格式与有限支持范围见 [原生 app](applications/prismspf/grain_growth_dream3d/README.md)。
声明使用 `.h`，实现使用 `.cpp`；应用入口沿 PRISMS-PF 使用 `.cc`。
本仓不包含 PRISMS-PF core、DREAM.3D、DAMASK、安装镜像或机器配置；产品通过自己的安装声明选择它们。

## 许可与来源

源码保留各文件 SPDX 与版权说明，采用 LGPL-2.1-or-later，全文见 [LICENSE](LICENSE)。
本应用依赖 [PRISMS-PF](https://github.com/PRISMS-Center/phaseField)，不替换上游 grainGrowth 示例。
上游 core 的许可与版权不因本仓分发而改变。见 [来源说明](THIRD_PARTY_NOTICES.md)。
