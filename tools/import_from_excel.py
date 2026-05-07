"""
从 TCMSP 手动导出的 Excel/CSV 文件导入数据到 tcmsp_builtin.json

使用方法：
1. 打开 https://www.tcmsp-e.com，搜索中药
2. 设置 OB≥30%, DL≥0.18 筛选
3. 点击页面上的导出/下载按钮，保存为 Excel 或 CSV
4. 运行：python3 tools/import_from_excel.py 黄连 /path/to/downloaded_file.xlsx
"""

import sys, json, pandas as pd
from pathlib import Path

BUILTIN_PATH = Path(__file__).parent.parent / "data" / "tcmsp_builtin.json"


def import_file(herb_name: str, file_path: str, ob_th: float = 30.0, dl_th: float = 0.18):
    path = Path(file_path)
    if not path.exists():
        print(f"❌ 文件不存在：{file_path}")
        return

    # 读取文件
    if path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    elif path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
    else:
        print("❌ 仅支持 .xlsx / .xls / .csv 格式")
        return

    print(f"原始数据：{len(df)} 行，列名：{df.columns.tolist()}")

    # 自动识别列名
    col_map = {}
    for c in df.columns:
        lo = c.lower().strip()
        if any(k in lo for k in ["mol_name", "molecule", "name", "成分"]):
            col_map["mol_name"] = c
        elif lo in ("ob", "ob (%)") or "oral" in lo:
            col_map["OB"] = c
        elif lo in ("dl", "drug") or "likeness" in lo:
            col_map["DL"] = c
        elif lo in ("mw", "mol_weight") or "weight" in lo:
            col_map["MW"] = c
        elif "smiles" in lo:
            col_map["SMILES"] = c

    if "mol_name" not in col_map:
        # 猜第一列是成分名
        col_map["mol_name"] = df.columns[0]
    if "OB" not in col_map or "DL" not in col_map:
        print("⚠️  无法自动识别 OB/DL 列，请手动检查列名：", df.columns.tolist())
        return

    # 转换并筛选
    df = df.rename(columns={v: k for k, v in col_map.items()})
    df["OB"] = pd.to_numeric(df["OB"], errors="coerce").fillna(0)
    df["DL"] = pd.to_numeric(df["DL"], errors="coerce").fillna(0)
    filtered = df[(df["OB"] >= ob_th) & (df["DL"] >= dl_th)]

    rows = []
    for _, row in filtered.iterrows():
        entry = {
            "mol_name": str(row.get("mol_name", "")).strip(),
            "OB":  round(float(row["OB"]), 2),
            "DL":  round(float(row["DL"]), 2),
        }
        if "MW" in row:
            entry["MW"] = round(float(row.get("MW") or 0), 2)
        if "SMILES" in row and str(row.get("SMILES", "")).strip():
            entry["SMILES"] = str(row["SMILES"]).strip()
        rows.append(entry)

    if not rows:
        print(f"⚠️  筛选后无数据（OB≥{ob_th}%, DL≥{dl_th}）")
        return

    # 写入 JSON
    if BUILTIN_PATH.exists():
        with open(BUILTIN_PATH, encoding="utf-8") as f:
            db = json.load(f)
    else:
        db = {}

    db[herb_name] = rows
    with open(BUILTIN_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, ensure_ascii=False, indent=2)

    print(f"✅ {herb_name}：导入 {len(rows)} 个活性成分")
    print(f"   已保存到 {BUILTIN_PATH}")
    print(f"\n下一步推送：")
    print(f"   git add data/tcmsp_builtin.json && git commit -m '补充{herb_name}内置数据' && git push")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法：python3 tools/import_from_excel.py <中药名> <文件路径>")
        print("示例：python3 tools/import_from_excel.py 黄连 ~/Downloads/huanglian.xlsx")
        sys.exit(1)

    herb = sys.argv[1]
    fpath = sys.argv[2]
    import_file(herb, fpath)
