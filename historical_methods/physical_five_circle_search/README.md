PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

# 显式物理五点共圆：8 条路径全部完成

本轮 **8 条新的物理五圈路径全部通过真实五点共圆判据，并越过至 s=1.02**，共 184 个接受节点、零校正失败。旧 e0093 + e0094 路线在原始 incircle 参数化下的事件后校正失败，已由此次稳定距离参数化解决；没有改进现有覆盖上界，完整精确上界仍为 `U=0.1155018241946`。

独立复验重新评分 **250 个保存坐标数组**，与原记录评分差全部为零。冻结 `event_forced_paths/` 的 702 文件和 `paired_event_search/` 的 655 文件全部保持原哈希。本目录是新的后续本地研究检查点，不更改此前 GitHub 首次归档树。

| 路径 | 事件对 | 接受节点采样峰值 R | 撤去全部附加约束后的事件后恢复 R |
|---|---|---:|---:|
| 00 | e0093 + e0094 | 0.11696087614379246 | 0.11684882799454933 |
| 01 | e0007 + e0072 | 0.1163859275528872 | 0.11637477030548703 |
| 02 | e0006 + e0002 | 0.11676527907073532 | 0.1167362542385903 |
| 03 | e0188 + e0235 | 0.11677028479842343 | 0.11674164328845031 |
| 04 | e0117 + e0065 | 0.11666630042327962 | 0.11656538295066544 |
| 05 | e0203 + e0050 | 0.11655836005239906 | 0.11655836005239896 |
| 06 | e0200 + e0038 | 0.1163518937397624 | 0.11631925034848167 |
| 07 | e0211 + e0212 | 0.11685999132724076 | 0.11681619136679153 |

本轮以 **203 个变量**直接表示几何：固定旋转片的 199 个中心自由度、平方覆盖半径 t、辅助圆心 q 的两个坐标、辅助平方圆半径 rho²。共同三点满足 `||ci-q||²=rho²`，两个 extra 点分别满足 `||cj-q||²=rho²+(1-s)Δj`，其中 Δj 是原共同三点最近空圆上的严格正平方距离 gap。s=1 时五点确实位于同一个辅助圆；s>1 时两个 extra 点进入原共同三点辅助圆，更新后的实际孔另由全盘 oracle 评估，不能宣称五点此时仍共圆。

防退化限制为共同三角形有向面积至少原值的 0.08，固定初始五点凸包循环顺序的每条边与其它三个顶点的 15 个有向面积至少原值的 0.05，`||q||≤0.999`、`0.03≤rho≤0.125`，并要求 95 个 nonunion 点位于辅助圆外。这些是此次有限搜索的几何限制，没有将限制内不可达解释为全局不可行。8 条路径接受节点的共同三角形面积比最小约 0.6533，循环边面积比最小约 0.3783，远高于所加下限；没有原始行列式由重合或共线造成的伪共圆。

筛查先要求两目标为真实盘内内部事件、共有三点、初始五点 Hull 恰为五顶点，初始三点空圆、两个 Δ 严格正及循环凸向可行。筛除旧 12 个 primitive-pair 路线（保留唯一合法 weak 五点案例）、并排除含 65 号盘的 cavity，避免重复另一个专门 65 号盘迁移任务。剩余 211 个候选只进行局部几何/Gram 筛查，取 weak 案例和 7 个低预测代价的新 cavity；只优化了表中 8 条路径。

初始响应使用完整 48 维 tangent 和平方半径目标 `G=ZᵀHZ`，以两个真实物理 gap 的归一化梯度构造 2×2 Gram，所有选定 Gram 正定、秩 2。控制列按 `extras` 的次序排列，事件 ID 对仅用于标识，不要求与该次序一致。包括全部 mixed 项的二阶响应同时校正原接触等式、两个 physical gap 等式及三个圆心/半径辅助等式。所有初始 5×203 圆等式 Jacobian 秩为 5；原 152 接触行加五个物理等式后的 157×203 Jacobian 满行秩 157。一阶物理响应残差不超过 `1.34e-15`，完整二阶物理响应残差不超过 `1.45e-15`。解析 Jacobian 的数值方向有限差检查绝对误差约为等式 `3.85e-10`、不等式 `2.13e-8`。

有限优化仍以 t 为目标，并动态更新完整实际孔的 `gk(C)≤t` 约束。辅助圆及其空圆约束只是几何控制，**不是覆盖半径 oracle**。每个保存节点用独立 `true_cover` 检查全部 100 个中心、4950 对的单位圆交点、所有边界反向驻点和盘内 Delaunay 顶点；全初始 DT、圆心穿界和边界切换库存也保留。所有 184 接受节点重新验证物理等式、非退化条件、circle Jacobian rank 5、q/中心/半径预算以及坐标哈希。

8 个 s=1 五点圆均在闭单位盘内、凸、最近空圆；五点平方距离散差在约 `10^-17～10^-16`。每条路径两事件均保存事件前后 `(+1,-1)` 的精确有理 incircle 符号，精确性仅针对保存的 binary64 坐标。所有路径都到达 s=1.02，每条 23 个接受节点。事件后先在原全孔问题中保留新侧，再降侧余量为零，最终移除所有附加圈、side 和 offset 条件，以原全 100 圆覆盖问题恢复和评分。

已知 `0.115516` 接触族的防重复检测同时要求 99 个 binding 中心的骨架对齐、active hole keys 匹配和独立 **全 100 中心** oracle 半径匹配。99 binding 指接触骨架，不表示 99 圆可以完成覆盖，也不表示 65 号盘可删。65 号盘始终参与实际孔、nonunion 空圆和 whole-D 评分；其非 binding 漂移只可能位于有界 rattler 区间，仍有 inactive cavity/slack 负债。本轮没有复制 65 号盘专向搬运搜索。

数学范围与后续限制：

- 本轮解决了稳定物理事件参数化及弱案例越过事件的问题，未找到更小 U。
- 表中峰值是保存的接受节点最大值，不能充当整段路径上确界证书。
- 五点共圆及空圆/inside 判据是数值验证；保存 binary64 点的精确符号不能认证代数事件零点或整段同伦。
- 局部 Gram 和 mixed 响应只适用于充分小控制变化、规则几何层，没有远事件最优性、严格局部极小或全球最优性证明。
- 物理圈半径、非退化和 empty 限制使此次路径更稳定，但也限制了可探索的几何；没有把数值失败或预算限制解释为全局下界。

文件：`input_checks.json`、`selected_cavities.json`、`screened_cavities.json`、`screen_rejections.json`、各 `predictorNN.npz` 保存筛查和 full mixed 局部数据；各 `pathNN_event1_event2/` 保存种子、轨迹、aux circle 状态、真实孔/事件库存、恢复和独立复验；`audit_report.json`、`candidate_hashes.json`、`checkpoint.json`、`MANIFEST_SHA256.json` 保存完成状态与哈希。`physical_paths.py` 为本轮源程序，`audit_physical.py` 只复验保存结果，不重跑优化。
