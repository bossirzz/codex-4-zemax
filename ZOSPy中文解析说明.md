# ZOSPy 中文解析适配

本机使用 ZOSPy 2.1.5。OpticStudio 中文报告能够正常读取，但部分高层分析封装使用英文专用语法：

- `SingleRayTrace` 期待 `Units`、`Coordinates`、`Wavelength` 和英文表格标题。
- `RayFan` 使用英文正则表达式匹配 `Tangential/Sagittal fan, field number ...`。
- `SystemData` 使用英文 Lark 语法；通用字段名称语法也不包含中文字符。

因此，之前报错发生在已读入报告之后的数据解析阶段。改文件编码不能消除这种字段语言差异。

`zospy_zh.py` 提供同名分析类，将已识别的简体中文标签转换成上游解析器要求的字段名。光扇仍直接读取原生数值数组。文件路径、题目、用户注释和光学数值不做全局翻译；缺失字段和未知格式仍报错。

从本目录运行脚本时用：

```python
from zospy_zh import SingleRayTrace, RayFan, SystemData

ray = SingleRayTrace(py=0.5).run(oss, oncomplete="Close")
fan = RayFan(number_of_rays=10).run(oss, oncomplete="Close")
system = SystemData().run(oss, oncomplete="Close")
```

项目的 `zospy_function_test.py` 已改用这三个适配类。直接调用原版 `zp.analyses` 的其他脚本仍使用原版解析器。Anaconda、当前虚拟环境中的 ZOSPy 安装及 OpticStudio 界面语言均未更改。

适配范围：本机 OpticStudio 24R1、ZOSPy 2.1.5，已保存的简体中文 DirectionCosines 单光线追迹、光扇和系统报告。其他追迹类型、不同版本报告、不同 Ray Aiming 状态以及其他分析的中文格式尚未验证。

验证文件：

- `zospy_zh_results/offline_report.json`：三项均通过，51 个单光线追迹数值单元与原始文本一致，82 个光扇数值与原生 CSV 一致，系统报告通过原版完整结果模型的必需字段校验。
- `zospy_zh_results/live_report.json`：重新连接 Extension 后，`SingleRayTrace`、`RayFan`、`SystemData` 三项现场验证均通过；测试前后镜头快照一致。结果保存在同目录的三个同名 JSON 文件中。没有保存镜头，测试结束已断开连接并保留 OpticStudio 窗口。

离线回归验证：

```powershell
& 'D:\workpalce\aizemaxp\zospy-python-env\.venv\Scripts\python.exe' -B '.\test_zospy_zh.py'
```

启用 OpticStudio Interactive Extension 并确认 API 授权可用后，可追加 `--live` 进行现场验证。
