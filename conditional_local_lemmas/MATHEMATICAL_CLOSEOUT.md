# Scoped local mathematics: public closeout map

The current universal result is .11323 < r_D(100) <= .1154469548855829;
see the sibling certified_bounds_11323 package. Exact global optimality remains
open. This is a curated exposition of earlier independently reviewed local
mathematics, not a new universal theorem or mathematical replay.

Private migration, Library identifiers and launcher locations are omitted.
The single latest canonical research tree supplies exact source/data, paired
reviews, counterexamples, corrections and explicitly cancelled drafts.

## 已审成果及精确作用域

| 成果 | 已接受结论 | 作用域及证据 |
|---|---|---|
| 实际三面核心 | fee>.025 | 同一实际共享apex，正向三面、circ≤R、covered unit anchors；完整源/反例/独审/两个不同隔离入口复现都在core包 |
| 三面更强及连续推广 | fee>.038；连续R≤794400/5378209 | actual signed-A；固定中央参数；独审源及精确端点守卫在extensions包 |
| 三面宽参数 | λ∈[.126,.154], t0±.014, R≤794400/5378209 ⇒ fee>.0205 | 匹配完整support；受限曲率移动cap，不是假设整个acute域凹；broad包保留失败算术史 |
| 旧弱两面路线 | 极弱fee及两个局部反例 | source-audit包的接受边界严格保留；不能当strongarc域反例或global反例 |
| 旧强两面 | τ=0 ⇒ .0223；generalτ ⇒ .00008 | stronggeometry包对应各自actual norm/arc条件，不自动扩成fullτ .0223 |
| 完整实际τ local | 固定R=.11322 ⇒ D>.0223 | SAME actual中间norm≤sqrt(1−R²)、两侧actual connected fullarc>.147；FULL_ACTUAL_TAU_PROOF及所有clarifications、完整独审 |
| 原Pro实际τ | 固定R=.11322 ⇒ D>.0225 | SAME actual中间norm²≤1−R²、两侧actual norm²≥.83；48exactguards/8whole-slabs历史ONE根复现JSONMATCH；原源不改 |
| 连续R local | R∈[.11322,.1135] ⇒ D>.0223 | 中间actual norm≤1、两侧original connected fullarc>.147；独审PASS，author/reviewer math0 |
| 连续R Pro | R∈[.11322,.116] ⇒ D>.0225 | 中间actual norm²≤1.0015、各outer actual norm²≥.83；8whole-slab解析radius转移、tau支付.228<.23；独审PASS，math0 |
| 固定R22参数 local | narrow±.00002 fee>.022284554，wide fee>.02202554 | local原actual cap+两侧arc条件；匹配新support和单budget，独审PASS |
| 固定R22参数 Pro | λ±.0003,t±.0002 fee>.0222855；balanced±.0002 fee>.0223455 | Pro原norm条件，单独重证.025 bootstrap；独审PASS |
| 新均值/弧和引理 | low D≤.0225 ⇒ M<H_R(D)；总弧≥.294至R=.11364或总弧≥.3018至R=.116 ⇒ D>.0225 | 中间同actual norm²≤1.0015；无outerfloor前提，两个不同actualouter各自connectedfullarc作Jensen；新paper独审PASS |
| 新R×参数组合 | Dλ,t>.0225−.6abs(Δλ)−.1725abs(Δt) | 上述同一actualfan弧和域，a≤.0003,b≤.0002；continuousR组合在新paper独审PASS |

新均值paper `q2_tau_mean_arc_envelope/PROOF.md` SHA256 `570696466a5e90d0d06e1738e41a4a94571ec7a8c45a66559fea27f2bd3cd663`；独审 `q2_tau_mean_arc_envelope_review/INDEPENDENT_REVIEW.md` SHA256 `7dea3c40dff003d166c9a26c39204d1f88f85dd65ea37c058d9bd1bcc6142390`。公式是 H_R(d)=1+R²−2R cos[.522+(100/99)(R−.11322)+.71sqrt(d)]，仅0≤d≤.0225。不等弧和闭阈值成立，实际同τ预算只付一次，不消耗中间HALFarc variance/Q/独立J。

`.147`单弧能推出outer norm²≥.83的临界R在(.113281,.113282)。这是旧保守接口瓶颈，不能误写成fee本身到此失效；新相关均值比较保留总弧.294至.11364。Pro/local skew两套约定不同：Pro alpha=3s/2+Z/2+eta, h=s/2−Z/2−eta；local同时反转d符号，公式不能混搭。

## 未审候选、未运行代码和保存状态

`.0228` cutoff仅未审 feasibility note：`q2_cutoff_headroom/UNREVIEWED_CANDIDATE.md`，SHA256 `f4551822a9fef0fb205302c3240f55d855e9814c39dc4941fe4ac43c7c3867b3`。旧slabs未覆盖新s/Delta/omega域；提出angular Lipschitz转移但未独审。特别旧signed tangent-gradient的7/40、23/20常数使用abs(Delta)≤.15，扩大到.151后不能照搬，改域sign仍待证。**没有D>.0228 theorem或全局新下界。** `.023`仅开始比较，没有结论。按收尾停止，不继续闭合。

`q2_uniform_fresh_guard/verify.py` SHA256 `d65a5094576b9851755720ca38ca8f2a1c06a47e9f66eab64e3923471b7868e7`是新常数代码草稿，**未静态代码独审、未执行、未编译、未导入**。新paper的有理比较已纸面独审；代码PASS没有发生。其GUARD_PLAN中的未来运行建议已被用户收尾指示撤销。本次连续推广/newmean/root/freshcode数学执行均0，旧Pro48和大链未重跑。

`q2_global_admission_interface`保留scope提取、投影覆盖规范化、最优配置存在性及明确限制。投影引理有child已完成纸面检查记录：投影外圆心到闭单位盘不会丢失旧同label已覆盖点/固定弧点集，但nearestowner/canonicalarc/结构fan/circ不能继承。反例原triangle (.98,.02),(.99,.01),(1.02,0) area1/10000,circ1/20，投影第三点后共线。覆盖半径在标号圆心supmetric下1-Lipschitz及紧致性给最优配置存在，不给数值。该child检查的范围/状态单独记录，不冒充新的整链独审。

## 失败路线与必须保留的纠错历史

| 路线/错误 | 失败原因与正确边界 | 应保留证据 |
|---|---|---|
| 把独立J fee套到真实A | J−A即真实tau，不能漏付或双花；显式37独立J障碍已审 | core包完整explicit37/attempt/review史 |
| 弱arc xi>.6推出大q2 fee | sin(epsilon)=1/12的实际连续/Delaunay/nearest局部点，Dactual=x/36<旧delta | q2_counterexample、ADDITIONAL_SCOPE、independent_audit、q2_compatibility |
| joint variance旧域 | cosgamma=17/25、tan(skew/2)=1/40，真实loss+中间HALFarc一次variance<delta；outerxi≈.592104，Gamma=1.65 | q2_joint_variance及独审；不能外推为fullarc>.147/Gamma1.85反例 |
| q3 broad首算 | 错误5eps/curvature=7/388；正确35/1944 | q3_parameter_broad/attempts/author_001与后成功receipt/traceback |
| τ0/mean/smallτ严格零乘积 | strict coefficient×zero不能仍strict；独立正margin供最终strict，另存clarifications | 三套ISSUE_HISTORY与所有ENDPOINT/ZERO_B/ZERO_SLACK补充 |
| 仅凭|d|<pi/2作skew身份 | 需先actual loss给|d|<.18及cos2eta>0，再付common-skew | SKEW_LOCALIZATION_CLARIFICATION及fulltau独审 |
| 跨R继承.522 baseline/.83 floor | baseline随R变，单弧normfloor随R降低；continuous derivation不可用点采样替代 | q2_radius_uniform_tau及geometrythresholds/current+coarsehistory |
| scalar factor2旧版 | 历史错误版本不能当现有conjugate；current proof/cert及独审才有效 | dependency research/quadratic_lambda14_conjugate/history_factor2_correction及current精确pins |
| 中间cap投影/反射未经fan重证 | covering投影成立但结构fan可能共线、circ/score改变 | q2_global_admission_interface及原sameactualcap scope |
| 启动/patch/xattr/materialize失败 | 非数学失败须与数学失败分开；不能据失败重做已成功数学 | 原attempt/receipt保留；library回执中array纠正、filepath材化纠正 |

这些失败史是防止重复走错的研究资料，应进入仓库的相关canonical history目录。重复制品、临时解压树和工具运行目录不需要进入仓库；本地原件保持不删。

## 复现入口与当前停止授权

已发布包的README、sourcepins、package manifest和历史execution receipts是身份/作用域入口。该原Pro数学源的公开相对路径入口 `research/pro_tau_local_replay_20261005/verify.py`运行会写同目录verification.json，历史根运行已在isolated-copy收到48/8 PASS且JSONBYTE_MATCH；目前不要运行。q3 core的两个历史isolated复现、本地Fraction端点、weak/joint反例入口均保留原运行源与路径限制，未为了搬运改checker绝对路径。不把静态ZIP解压PASS说成数学relocated replay。

原包的封装工具与delivery/acquisition记录属于历史材料，本公开子集未复制；当前公开文件由PUBLIC_COPY_RECEIPT和MANIFEST_SHA256绑定。新freshguard及候选不在已审复现成功清单。按当前用户指示，所有数学复现都暂停，不增加批次。

