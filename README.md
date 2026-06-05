# Robust PCA Seismic Denoising

这个仓库包含两份 MATLAB 单文件脚本，用于复现和扩展 Cheng, Chen & Sacchi (2015) 启发的 f-x-y Robust PCA / Principal Component Pursuit 地震随机强噪声压制流程。

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
| `src/rpca_tapered_reproduction.m` | 较简洁的 tapered f-x-y RPCA/PCP 复现实验，包含合成数据、平滑频率 taper、PCA 基线、RPCA 和参数扫描。 |
| `src/rpca_engineering_workflow.m` | 工程化入口，支持滑动时间窗 overlap-add、自适应频率相关 lambda、自适应 rank 投影、MAT 输入和可选 SEG-Y I/O 钩子。 |

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

- 合成数据实验可以计算 `Q_full` 和 `Q_tapered_ref` 指标；真实数据没有 clean reference 时会跳过 Q 指标。
- `rpca_engineering_workflow.m` 更适合真实数据流程；`rpca_tapered_reproduction.m` 更适合算法对照和参数理解。
- 单文件脚本保留了局部函数，方便直接运行和迁移。如果后续项目继续扩大，可以再拆分为 `functions/` 或 MATLAB package。
