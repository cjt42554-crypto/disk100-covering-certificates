# disk100严格局部极小性补充证书

已精确证明冻结的 `improved_algebraic_1154469` 100圆覆盖构型为**同时旋转商空间中的严格局部极小值**。半径约为 `0.11544695488558220163979840676978324705424`。这不是全局最优证书。

独立标准库接受：

```sh
python3 -S -O verify_sosc.py
python3 -S -O test_negative.py
```

从任意工作目录调用这两个脚本也可，路径以脚本本身为准。接受程序重新计算932变量KKT根箱、100圆全闭单位圆盘覆盖、孔的隐函数资格，以及439×493几何约束Jacobian的54维切空间Hessian严格正定性。数值提议文件和数值特征值均不作为接受依据。

`verification.json` 保存冻结时的精确有理界，重放接受将新报告写为 `replay_verification.json` 并打印stdout。经dyadic W和真零空间基Z变换的54×54 Hessian，其最小特征值的严格下界约为 `0.9999808203135334`；对应总误差上界约 `0.0000191796864665582<1`。这里的下界不属于原Hessian或正交切空间矩阵。标准库完整接受约11秒，SOSC整数步骤约1秒。

`negative_tests.json` 保存8项损坏拒绝。`SOURCE_BINDING.json` 绑定冻结源manifest SHA `4bc4d4647c5b5e8e6e7bb01d1ce10f50fae1823139ce0333efffa15f329472e2` 和实际复制的8个依赖SHA。`source/`是独立副本，上游冻结目录未被修改。接受器只在自身目录写再生报告。

`LOCAL_PROOF.md` 给出Neumann零空间证书、严格正定接受不等式、146个真实接触不等式与439条等式KKT的桥、辅助孔与真正全闭盘覆盖目标的关系，以及旋转规范的严格局部解释。`numeric_diagnose.py` 和 `propose.py` 是已从本公开子集省略的历史NumPy/SciPy数值提议工具；标准库重放不依赖它们。

本公开副本只作静态整理与metadata隐私清理，未新增数学重放。接受源码及精确witness逐字节未变；原执行证据保留为历史身份。
