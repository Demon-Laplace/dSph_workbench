# 无银河系引力的理论来流风洞

使用项目登记的投放初态、MW1 预期剖面来流和冻结风洞求解器，保留内部旋转、自引力、气体物理和共同参考系变换。维护实现位于科研项目 Code/simulation/wind_tunnel/；这里是同步的服务器运行源码。

run_server_wind.py stage 只打包已有资料，不启动本机模拟；显式传入 --project、--source-runtime、--output、--manifest。运行包复用项目 Code/analysis/model_comparisons/wind_tunnel_response.py，输出仅含该模拟自身的时域测量及图件。

上传后在独立服务器目录使用 wind_only_astrothin.slurm：astro_thin 单节点24 MPI、每进程单线程、Intel MPI shm与核绑定、72000 MiB、99小时。编译和积分均在计算节点，运行至少140秒且时间推进后生成startup.json；这不是完整科学验证。默认从同一个初态启动；原生检查点只在原兼容构建和24 MPI布局中续跑，续跑时以脚本参数1代替默认2。当前节点没有GPU资源。

TimeMax从首步使用完整目标值，短段停止使用原生output/stop或墙钟限制。实验输出不覆盖旧轮；仅按本次用户授权提交或续跑。
