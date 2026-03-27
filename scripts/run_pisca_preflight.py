from __future__ import annotations

import json
import os
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PISCA_DIR = ROOT / "PISCA"


def _exists(p: Path) -> bool:
    return p.exists()


def _check_beast_root() -> dict[str, object]:
    beast_root_env = os.environ.get("BEAST_ROOT", "").strip()
    candidates: list[Path] = []
    if beast_root_env:
        candidates.append(Path(beast_root_env))
    # Common locations users may place BEAST 1.8.x
    candidates += [
        Path(r"D:\tools\BEASTv1.8.4"),
        ROOT / "tools" / "BEAST v1.8.4",
        Path(r"D:\BEASTv1.8.4"),
        Path(r"C:\tools\BEASTv1.8.4"),
        Path.home() / "BEASTv1.8.4",
    ]

    found: Path | None = None
    for c in candidates:
        has_bin = (c / "bin" / "beast").exists() or (c / "bin" / "beast.cmd").exists()
        if (c / "lib" / "beast.jar").exists() and has_bin:
            found = c
            break

    return {
        "beast_root_env": beast_root_env or None,
        "beast_root_found": str(found) if found else None,
        "beast_jar_exists": bool(found and (found / "lib" / "beast.jar").exists()),
        "beast_bin_exists": bool(
            found and ((found / "bin" / "beast").exists() or (found / "bin" / "beast.cmd").exists())
        ),
    }


def _check_pisca_plugin(beast_root: str | None) -> dict[str, object]:
    dist_jar = PISCA_DIR / "dist" / "PISCA.PISCALoader.jar"
    ant_cmd = shutil.which("ant")
    if not ant_cmd:
        bundled_ant = ROOT / "tools" / "apache-ant-1.10.15" / "bin" / "ant.bat"
        if bundled_ant.exists():
            ant_cmd = str(bundled_ant)
    java_cmd = shutil.which("java")
    jar_cmd = shutil.which("jar")

    installed_plugin = None
    if beast_root:
        p = Path(beast_root) / "lib" / "plugins" / "PISCA.PISCALoader.jar"
        if p.exists():
            installed_plugin = str(p)

    return {
        "pisca_source_dir": str(PISCA_DIR),
        "pisca_dist_jar_exists": dist_jar.exists(),
        "pisca_dist_jar": str(dist_jar),
        "pisca_installed_plugin": installed_plugin,
        "java_in_path": bool(java_cmd),
        "jar_in_path": bool(jar_cmd),
        "ant_in_path": bool(ant_cmd),
    }


def _check_required_data() -> dict[str, object]:
    # Strictly matching the source Rmd path used by authors for Fig4AB scatter/methylation source
    fig4ab_source = ROOT / "Revision" / "Results" / "Longitudinal_meth_data" / "Data_source_methylation_Fig.4_SCLL12-SCLL19.xlsx"
    moesm8 = ROOT / "evoflux" / "41586_2025_9374_MOESM8_ESM.xlsx"

    return {
        "fig4ab_source_expected": str(fig4ab_source),
        "fig4ab_source_exists": fig4ab_source.exists(),
        "moesm8_path": str(moesm8),
        "moesm8_exists": moesm8.exists(),
    }


def _required_actions(report: dict[str, object]) -> list[str]:
    actions: list[str] = []
    beast = report["beast"]  # type: ignore[index]
    plugin = report["pisca_plugin"]  # type: ignore[index]
    data = report["data"]  # type: ignore[index]

    if not beast["beast_root_found"]:
        actions.append("安装 BEAST 1.8.4（必须 1.8.x），并设置环境变量 BEAST_ROOT 指向其根目录。")
    if not plugin["pisca_installed_plugin"]:
        if plugin["pisca_dist_jar_exists"]:
            actions.append("执行 PISCA/install.sh <BEAST_ROOT> 安装已编译插件。")
        else:
            actions.append("当前 PISCA 为源码仓（无 dist jar）；请在 PISCA 下配置 beast_sdk.properties 后运行 `ant test-install` 编译并安装插件。")
    if not data["fig4ab_source_exists"]:
        actions.append("补充作者 Fig4AB 原始输入：Revision/Results/Longitudinal_meth_data/Data_source_methylation_Fig.4_SCLL12-SCLL19.xlsx。")
        if data["moesm8_exists"]:
            actions.append("若无法提供上述文件，可先用 MOESM8 构建替代输入（可跑流程，但与论文 Data_source 口径有差异）。")
    return actions


def main() -> None:
    beast = _check_beast_root()
    plugin = _check_pisca_plugin(beast.get("beast_root_found"))  # type: ignore[arg-type]
    data = _check_required_data()

    report = {
        "root": str(ROOT),
        "beast": beast,
        "pisca_plugin": plugin,
        "data": data,
    }
    report["required_actions"] = _required_actions(report)

    out_json = ROOT / "figures" / "reproduced" / "pisca_preflight_report.json"
    out_md = ROOT / "figures" / "reproduced" / "pisca_preflight_report.md"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines: list[str] = []
    lines.append("# PISCA 复现预检查报告")
    lines.append("")
    lines.append(f"- 工作目录：`{ROOT}`")
    lines.append("")
    lines.append("## BEAST")
    lines.append(f"- `BEAST_ROOT` 环境变量：`{beast['beast_root_env']}`")
    lines.append(f"- 发现 BEAST 根目录：`{beast['beast_root_found']}`")
    lines.append("")
    lines.append("## PISCA 插件")
    lines.append(f"- 源码目录：`{plugin['pisca_source_dir']}`")
    lines.append(f"- dist jar 存在：`{plugin['pisca_dist_jar_exists']}`")
    lines.append(f"- 已安装插件：`{plugin['pisca_installed_plugin']}`")
    lines.append(f"- java in PATH：`{plugin['java_in_path']}`")
    lines.append(f"- ant in PATH：`{plugin['ant_in_path']}`")
    lines.append("")
    lines.append("## 数据")
    lines.append(f"- Fig4AB Data_source 文件：`{data['fig4ab_source_expected']}`")
    lines.append(f"- 是否存在：`{data['fig4ab_source_exists']}`")
    lines.append(f"- MOESM8 文件：`{data['moesm8_path']}`")
    lines.append(f"- 是否存在：`{data['moesm8_exists']}`")
    lines.append("")
    lines.append("## 下一步必做")
    if report["required_actions"]:
        for i, a in enumerate(report["required_actions"], start=1):
            lines.append(f"{i}. {a}")
    else:
        lines.append("1. 预检查已通过，可进入 XML 构建与 BEAST 运行。")
    lines.append("")
    out_md.write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"\nWROTE {out_json}")
    print(f"WROTE {out_md}")


if __name__ == "__main__":
    main()

