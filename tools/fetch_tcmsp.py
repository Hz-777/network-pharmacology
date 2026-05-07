"""
本地运行脚本：从 TCMSP 批量抓取中药成分数据，保存到 data/tcmsp_builtin.json
用法：python tools/fetch_tcmsp.py
"""

import sys, json, time, requests
from pathlib import Path
from bs4 import BeautifulSoup

sys.path.insert(0, str(Path(__file__).parent.parent))

BUILTIN_PATH = Path(__file__).parent.parent / "data" / "tcmsp_builtin.json"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
}

# ── 在这里添加你想要抓取的中药名称 ──────────────────────────────────────────
HERBS_TO_FETCH = [
    "黄连", "黄芩", "甘草", "人参", "黄芪", "丹参", "当归", "川芎",
    "柴胡", "茯苓", "白术", "陈皮", "大黄", "麻黄", "桂枝", "附子",
    "半夏", "金银花", "连翘", "赤芍", "生地黄", "枸杞子", "菊花",
    "薄荷", "赤石脂", "粳米", "鳖甲", "干姜", "肉桂", "防风", "白芍",
    # 在这里继续添加 ↓
    "五味子", "酸枣仁", "远志", "石菖蒲", "龙骨", "牡蛎", "磁石",
    "天麻", "钩藤", "石决明", "珍珠母", "羚羊角",
    "牛黄", "麝香", "冰片", "苏合香",
    "三七", "红花", "桃仁", "益母草", "泽兰",
    "木通", "车前子", "泽泻", "茵陈", "金钱草",
    "杏仁", "川贝母", "浙贝母", "桔梗", "枇杷叶",
    "砂仁", "豆蔻", "藿香", "佩兰", "苍术",
    "葛根", "升麻", "蝉蜕", "牛蒡子",
    "板蓝根", "蒲公英", "紫花地丁", "鱼腥草",
    "熟地黄", "山药", "山茱萸", "泽泻", "牡丹皮",
]

OB_THRESHOLD = 30.0
DL_THRESHOLD = 0.18


def fetch_one(herb: str, session: requests.Session) -> list:
    """Fetch compounds for one herb from TCMSP."""
    for base in ["https://www.tcmsp-e.com", "https://old.tcmsp-e.com"]:
        try:
            # Try JSON API
            resp = session.post(
                f"{base}/api/getmolecules.php",
                data={"herb_cn_name": herb, "token": ""},
                timeout=20,
            )
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list) and data:
                    return _parse_json(data)
        except Exception:
            pass

        try:
            # Fallback: HTML scraping
            resp = session.get(
                f"{base}/tcmsp.php",
                params={"qr": herb, "qsr": "herb_cn_name"},
                timeout=20,
            )
            rows = _parse_html(resp.text)
            if rows:
                return rows
        except Exception:
            pass

    return []


def _parse_json(data: list) -> list:
    results = []
    for m in data:
        try:
            ob = float(m.get("ob") or m.get("OB") or 0)
            dl = float(m.get("dl") or m.get("DL") or 0)
            if ob >= OB_THRESHOLD and dl >= DL_THRESHOLD:
                results.append({
                    "mol_name": m.get("mol_name") or m.get("name") or "",
                    "OB":  round(ob, 2),
                    "DL":  round(dl, 2),
                    "MW":  round(float(m.get("mw") or m.get("MW") or 0), 2),
                    "SMILES": m.get("smiles") or m.get("SMILES") or "",
                })
        except Exception:
            continue
    return results


def _parse_html(html: str) -> list:
    soup = BeautifulSoup(html, "lxml")
    for table in soup.find_all("table"):
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        if not any(k in " ".join(headers) for k in ["Molecule", "OB", "DL"]):
            continue
        col = {h: i for i, h in enumerate(headers)}
        results = []
        for tr in table.find_all("tr")[1:]:
            tds = [td.get_text(strip=True) for td in tr.find_all("td")]
            if not tds:
                continue
            try:
                ob = float(tds[col.get("OB (%)", col.get("OB", 999))])
                dl = float(tds[col.get("DL",  999)])
                if ob >= OB_THRESHOLD and dl >= DL_THRESHOLD:
                    results.append({
                        "mol_name": tds[col.get("Molecule Name", 0)],
                        "OB": round(ob, 2),
                        "DL": round(dl, 2),
                    })
            except Exception:
                continue
        if results:
            return results
    return []


def main():
    # Load existing data
    if BUILTIN_PATH.exists():
        with open(BUILTIN_PATH, encoding="utf-8") as f:
            db = json.load(f)
        print(f"已有内置数据：{len(db)} 味中药")
    else:
        db = {}

    session = requests.Session()
    session.headers.update(HEADERS)

    new_count = 0
    for i, herb in enumerate(HERBS_TO_FETCH):
        if herb in db:
            print(f"  [{i+1}/{len(HERBS_TO_FETCH)}] {herb} — 已存在，跳过")
            continue

        print(f"  [{i+1}/{len(HERBS_TO_FETCH)}] 正在抓取：{herb} ...", end=" ", flush=True)
        rows = fetch_one(herb, session)

        if rows:
            db[herb] = rows
            new_count += 1
            print(f"✅ {len(rows)} 个活性成分")
        else:
            print("⚠️  未获取到数据")

        time.sleep(1.5)  # 礼貌性延迟，避免被封

    # Save
    with open(BUILTIN_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, ensure_ascii=False, indent=2)

    print(f"\n完成！新增 {new_count} 味中药，共 {len(db)} 味。")
    print(f"文件已保存到：{BUILTIN_PATH}")
    print("\n下一步：git add data/tcmsp_builtin.json && git push")


if __name__ == "__main__":
    main()
