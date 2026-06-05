# Robust PCA Seismic Denoising

这个仓库包含两份 MATLAB 单文件脚本，用于复现和扩展 Cheng, Chen & Sacchi (2015) 提出的 f-x-y 域 Robust PCA / Principal Component Pursuit 思路，目标是压制地震数据中的 erratic noise，也就是空间上离群、幅值可能很强、统计分布明显非高斯的随机强噪声。

## 论文来源

本仓库主要参考的原文是：

> Cheng, J., Chen, K., & Sacchi, M. D. (2015). **Application of Robust Principal Component Analysis (RPCA) to suppress erratic noise in seismic records**. *SEG Technical Program Expanded Abstracts 2015*, 4646-4651.

- DOI: <https://doi.org/10.1190/segam2015-5869427.1>
- PDF: <https://saig.physics.ualberta.ca/lib/exe/fetch.php?media=publications%3Aconfpro%3Aseg2015_cheng_2.pdf>

## 背景与目标

常规地震随机噪声压制方法通常假设噪声比较接近高斯分布，或者能在某个变换域中被平滑地削弱。但 field data 中经常会出现 erratic noise，例如坏道、脉冲干扰、采集系统异常、局部强能量污染等。这类噪声的特点是幅值大、位置稀疏、分布不规则，因此普通最小二乘类方法和固定 rank 的 f-x-y PCA/eigenimage 滤波容易受到离群值影响。

Cheng, Chen & Sacchi 的核心想法是：在固定频率切片上，有效地震信号由于空间相干性，往往可以用低秩矩阵近似；erratic noise 在空间上只污染部分位置，适合用稀疏矩阵表示。因此可以把每个频率切片分解为：

```text
D(f) = L(f) + S(f)
```

其中：

- `D(f)` 是某个频率下的输入 f-x-y 切片；
- `L(f)` 是低秩部分，对应相干地震有效信号；
- `S(f)` 是稀疏部分，对应 erratic noise 或强离群污染。

## 算法原理

### f-x-y 频率切片

输入数据一般写成 `[nt, nx, ny]`。代码先沿时间轴做 FFT，把数据从 t-x-y 域变到 f-x-y 域。对每一个正频率 `f`，取出一个空间矩阵：

```text
D_f = D(f, :, :)
```

二维数据 `[nt, nx]` 会被当成 `[nt, nx, 1]` 处理；三维数据则直接形成 `nx` by `ny` 的频率切片。

### PCA/eigenimage 基线

传统 f-x-y eigenimage filtering 对每个频率切片做 SVD：

```text
D_f = U Sigma V*
```

然后只保留前 `r` 个奇异值：

```text
L_f = U(:,1:r) Sigma(1:r,1:r) V(:,1:r)*
```

这种方法对较平稳的随机噪声有效，但当某些道或局部位置被强离群噪声污染时，SVD 的主成分也会被拉偏。因此本仓库保留 PCA 作为 baseline，用来和 RPCA 结果对比。

### RPCA / PCP 分解

RPCA 使用 Principal Component Pursuit 的凸优化形式：

```text
minimize    ||L||_* + lambda ||S||_1
subject to  D = L + S
```

这里：

- `||L||_*` 是 nuclear norm，即奇异值之和，用来鼓励 `L` 低秩；
- `||S||_1` 是逐元素 L1 范数，用来鼓励 `S` 稀疏；
- `lambda` 控制多少能量被分到稀疏噪声项中。

代码中使用 inexact augmented Lagrange multiplier / IALM 思路迭代求解。每次迭代主要包含两个阈值操作：

- 对奇异值做 soft-thresholding，得到低秩矩阵 `L`；
- 对矩阵元素做 soft-thresholding，得到稀疏矩阵 `S`。

求得每个频率切片的 `L_f` 后，代码按共轭对称关系补齐负频率，再做 inverse FFT 回到时间域。最终：

- `rpca_denoised` 是低秩重建结果；
- `rpca_sparse_noise` 是 RPCA 分离出的稀疏噪声；
- `removed_noise = noisy - rpca_denoised` 是工程上更直观的去除量。

## 本仓库的做法

原文强调在 f-x-y 域用 RPCA 把低秩有效信号和稀疏 erratic noise 分开。本仓库在这个主线基础上做了几个工程化处理，目的是让脚本更稳、更容易用于 synthetic test 和真实数据试跑。

### 1. 平滑频率 taper

代码没有使用硬的频率截断，而是用 raised-cosine spectral taper：

```matlab
freq_taper.flow_stop  = 0;
freq_taper.flow_pass  = 2;
freq_taper.fhigh_pass = 45;
freq_taper.fhigh_stop = 75;
```

`flow_pass` 到 `fhigh_pass` 是完整保留频带，低频和高频边界使用余弦过渡。这样做可以减少硬截频带来的 Gibbs/ringing artifact。

### 2. 时间轴补零

FFT 前使用 `pad_factor` 对时间轴补零，使频率采样更平滑，降低 inverse FFT 后的循环泄漏和边界振铃。

### 3. RPCA 后 rank projection

纯 PCP 已经能把 `D = L + S` 分开，但真实地震数据中还会有较弱的非相干噪声残留。本仓库在 RPCA 之后可选地对 `L` 再做一次 rank projection：

```matlab
opts.use_rank_projection = true;
opts.rank_mode = 'adaptive';
```

这不是原始 PCP 的必要步骤，而是一个实用增强：先用 RPCA 去掉强离群噪声，再用低秩投影压制残余随机能量。

### 4. 滑动时间窗 overlap-add

`rpca_engineering_workflow.m` 支持滑动时间窗：

```matlab
USE_SLIDING_WINDOW = true;
win_length_s = 0.50;
win_overlap  = 0.50;
```

每个时间窗内单独做 f-x-y RPCA，然后用正弦窗权重 overlap-add 拼回完整时间轴。这样可以缓解长记录中非平稳事件、频谱变化和局部噪声差异带来的问题。

### 5. 自适应 lambda 和 rank

工程化脚本支持频率相关的 `lambda(f)` 和 `rank(f)`：

```matlab
opts.lambda_adaptive = true;
opts.lambda_freq_boost = 0.40;
opts.rank_mode = 'adaptive';
opts.rank_min = 2;
opts.rank_max = 7;
```

直观上，低频有效信号通常更强、更相干，高频更容易受噪声影响。代码允许不同频率使用不同稀疏惩罚和 rank cap，以减少“一组参数管所有频率”的僵硬感。

### 6. 合成数据参数扫描

当 `INPUT_MODE = 'synthetic'` 时，脚本知道 clean reference，可以计算：

```text
Q = 10 log10(||clean||_F^2 / ||estimate - clean||_F^2)
```

并自动扫描 `lambda` 和 rank 参数。真实数据没有 clean reference 时，Q 指标会跳过，主要依赖 wiggle、image、removed noise 和 rank/iteration QC。

## 仓库结构

```text
.
├── src/
│   ├── rpca_tapered_reproduction.m
│   └── rpca_engineering_workflow.m
├── .gitignore
└── README.md
```

## 脚本说明

| 脚本 | 用途 |
| --- | --- |
| `src/rpca_tapered_reproduction.m` | 较简洁的 tapered f-x-y RPCA/PCP 复现实验，包含合成数据、平滑频率 taper、PCA 基线、RPCA 和参数扫描。适合理解算法主线。 |
| `src/rpca_engineering_workflow.m` | 工程化入口，支持滑动时间窗 overlap-add、自适应频率相关 lambda、自适应 rank 投影、MAT 输入和可选 SEG-Y I/O 钩子。适合真实数据试跑。 |

## 环境要求

- MATLAB，建议 R2018b 或更新版本。
- 不需要额外工具箱即可运行合成数据示例。
- SEG-Y 输入/输出依赖当前 MATLAB 环境中可用的 `segyread` / `segywrite`。

## 快速开始

在 MATLAB 中把当前目录切到仓库根目录，然后运行工程化示例：

```matlab
cd('path/to/Robust PCA')
run('src/rpca_engineering_workflow.m')
```

如果只想跑较轻量的复现实验：

```matlab
run('src/rpca_tapered_reproduction.m')
```

默认配置会生成合成 3D 地震数据，加入高斯噪声和 erratic noise，然后比较以下结果：

- smooth spectral taper reconstruction
- f-x-y PCA/eigenimage filtering
- enhanced f-x-y RPCA/PCP denoising
- RPCA removed noise / sparse noise QC

## 使用自己的 MAT 数据

工程化脚本中修改这些参数：

```matlab
INPUT_MODE = 'mat';
USER_DATA_MAT = 'your_data.mat';
DATA_VAR_NAME = 'data';
dt = 0.002;
```

`DATA_VAR_NAME` 对应的数组应为：

- `[nt, nx]`：二维道集
- `[nt, nx, ny]`：三维数据体或多个 y slice

脚本会自动把二维输入 reshape 为 `[nt, nx, 1]`。

## 关键参数

- `freq_taper`：平滑频率 taper，避免硬截频造成的 Gibbs/ringing artifact。
- `pad_factor`：FFT 时间轴补零倍数。
- `USE_SLIDING_WINDOW`、`win_length_s`、`win_overlap`：滑动时间窗处理参数。
- `opts.lambda_scale`：RPCA/PCP 稀疏项权重基准。
- `opts.lambda_adaptive`：是否启用频率相关 lambda。
- `opts.use_rank_projection`、`opts.rank_mode`：RPCA 后的 rank 投影策略。
- `DO_SYNTHETIC_SWEEP`：合成数据上的参数扫描开关；真实数据建议关闭。

## 输出

默认会保存 MATLAB 结果文件：

- `rpca_tapered_reproduction_results.mat`
- `rpca_engineering_workflow_results.mat`

这些结果文件、SEG-Y 文件和大型 MAT 数据默认不会进入 Git 版本控制。需要发布示例数据时，建议单独准备小型匿名样例，并在 `.gitignore` 中显式放行。

## 注意事项

- 本仓库是研究复现和工程试验脚本，不是对原文结果的逐图逐数值完全复刻。
- 合成数据实验可以计算 `Q_full` 和 `Q_tapered_ref` 指标；真实数据没有 clean reference 时会跳过 Q 指标。
- `rpca_engineering_workflow.m` 更适合真实数据流程；`rpca_tapered_reproduction.m` 更适合算法对照和参数理解。
- 单文件脚本保留了局部函数，方便直接运行和迁移。如果后续项目继续扩大，可以再拆分为 `functions/` 或 MATLAB package。
