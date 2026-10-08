# DREAM.3D 组织输入的晶粒长大 app

在原版 PRISMS-PF v4.0.1 core 上运行普通多晶粒 Allen–Cahn 长大；`order_parameter_count`
是运行参数，不需要按 count 重新编译。本切片为二维/三维 Q1、单 MPI rank、固定矩形网格和自然边界。

```sh
cmake -S applications/prismspf/grain_growth_dream3d -B build -DCMAKE_BUILD_TYPE=Release -DPRISMS-PF_DIR=<core>/lib/cmake/PRISMS-PF
cmake --build build --parallel 2
cmake --install build --prefix <install>
```

在输入目录运行 `<install>/bin/grain_growth_dream3d`。产品 adapter 将 DREAM.3D 或 EBSD 组织
转成 `feature_ids.txt`，头为 `prisms_feature_grid_v1 dimension nx ny nz grain_count`，随后是
按 X 最快、Z 最慢排列的稠密正整数晶粒编号。PRM 保存网格、count、方程常数及固定初始化参数。
原始取向由产品保存；通道编号不代表晶体取向或物理晶粒编号。

`modules/grain_structure/upstream/` 的六个头文件/实现文件来自 PRISMS-PF
[v2.4.2 固定提交](https://github.com/PRISMS-Center/phaseField/tree/150d85f74994eca7ea0d756d164e28eccef97b12)，
保持原算法：连通区域识别、包围球分配、通道重映射。旧版只提供这些方法，不作为构建核心。
`native.cpp` 负责 v4 生命周期、离散组织采样、扩散平滑及逐帧身份输出；不需要启用 VTK reader。
初始晶粒几何直接取输入单元 FeatureId，避免先转节点标签再取众数造成小晶粒缺失；
后续连通识别仍使用上游方法。该修正不保证小晶粒在平滑后的连续场中保有优势区域。
普通长大势阱来自 `modules/grain_growth/model.*`。这些适配不等于与旧版完整算例数值等价。

连续场输出为 VTU；同帧 `*-grains-*.json` 保存当前连通区域对应的晶粒编号及通道映射。
当 count 不小于输入晶粒数时，每个晶粒独占固定通道，不调用重分配或最近中心识别；
同一晶粒的不连通部分保留其场身份，输出按单元中心优势场生成晶粒标签和几何观察。
复用通道时身份传递仍沿用上游最近中心方法，并检查通道一致性与重复身份；歧义不会被当作可用的 DAMASK 输入。
复用通道下同一输入晶粒断开成多个区域的持续身份尚未验收。分配成功只证明当前配置可用；
贪心包围球分配失败也不能证明该 count 在数学上不可能。当前真实整域 count=9 仍有冲突。
不暗改 count、删除小晶粒或合并晶粒来凑足通道。

count 是用户选择的实际正整数通道数，不限于 6/9，也不静默截成某个上限；
当前没有跨网格、跨机器通用的最大可运行值。较少通道可能不能完成当前几何分配，较多通道增加内存和计算量。
初始平滑及有限元采样可能使小晶粒失去优势单元；完整身份表不等于所有输入晶粒在初态仍有优势区域。

原生与软件证据由产品验收记录引用；短时规则网格案例不构成材料标定、网格独立性或长时间演化验证。
本 app 及适配采用 LGPL-2.1-or-later，保留 PRISMS-PF 上游许可；不包含储能、回复或再结晶出生。
