"""
用 requests 直接从 TCMSP 批量抓取中药成分数据
比 Selenium 更快更稳定，无需浏览器

用法：python3 tools/fetch_tcmsp_requests.py
"""

import json, re, time
import requests
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

BUILTIN_PATH = Path(__file__).parent.parent / "data" / "tcmsp_builtin.json"
BASE_URL     = "https://www.tcmsp-e.com"
OB_THRESHOLD = 30.0
DL_THRESHOLD = 0.18

# Chinese → English/Latin name fallbacks for herbs that fail Chinese search
CN_TO_EN_FALLBACK = {
    "天麻": "Gastrodia Rhizoma",
    "远志": "Radix Polygalae",
    "益母草": "Herba Leonuri",
    "车前子": "Semen Plantaginis",
    "泽泻": "Rhizoma Alismatis",
    "茵陈": "Herba Artemisiae Scopariae",
    "金钱草": "Herba Lysimachiae",
    "杏仁": "Semen Armeniacae Amarum",
    "桔梗": "Radix Platycodonis",
    "砂仁": "Fructus Amomi",
    "藿香": "Herba Agastaches",
    "苍术": "Rhizoma Atractylodis",
    "葛根": "Radix Puerariae Lobatae",
    "升麻": "Rhizoma Cimicifugae",
    "板蓝根": "Radix Isatidis",
    "蒲公英": "Herba Taraxaci",
    "鱼腥草": "Herba Houttuyniae",
    "山药": "Rhizoma Dioscoreae",
    "山茱萸": "Fructus Corni",
    "牡丹皮": "Cortex Moutan",
}

# Herbs to fetch
HERBS = [
    "黄连", "黄芩", "甘草", "人参", "黄芪", "丹参", "当归", "川芎",
    "柴胡", "茯苓", "白术", "陈皮", "大黄", "麻黄", "桂枝", "附子",
    "半夏", "金银花", "连翘", "赤芍", "生地黄", "熟地黄", "枸杞子",
    "菊花", "薄荷", "赤石脂", "粳米", "鳖甲", "干姜", "肉桂",
    "防风", "白芍", "五味子", "酸枣仁", "远志", "石菖蒲",
    "天麻", "钩藤", "三七", "红花", "桃仁", "益母草",
    "车前子", "泽泻", "茵陈", "金钱草", "杏仁", "桔梗",
    "砂仁", "藿香", "苍术", "葛根", "升麻", "板蓝根",
    "蒲公英", "鱼腥草", "山药", "山茱萸", "牡丹皮",
]


def _make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                      "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "zh-CN,zh;q=0.9,en;q=0.8",
        "Referer": f"{BASE_URL}/tcmsp.php",
    })
    return s


def _get_token(session: requests.Session) -> str:
    r = session.get(f"{BASE_URL}/tcmsp.php", timeout=20)
    tokens = re.findall(r"value=['\"]([a-f0-9]{32})['\"]", r.text)
    return tokens[0] if tokens else ""


def _search_herbs(session, token: str, query: str) -> list:
    """Search for herbs by name, return list of {herb_cn_name, herb_en_name, herb_pinyin}."""
    r = session.get(
        f"{BASE_URL}/tcmspsearch.php",
        params={"qs": "herb_all_name", "q": query, "token": token},
        timeout=20,
    )
    m = re.search(r"data:\s*(\[(?:[^\[\]]*|\[[^\[\]]*\])*\])", r.text, re.DOTALL)
    if not m or len(m.group(1).strip()) <= 2:
        return []
    try:
        return json.loads(m.group(1))
    except Exception:
        return []


def _get_compounds(session, token: str, herb_en_name: str) -> list:
    """Get all compounds from herb detail page."""
    r = session.get(
        f"{BASE_URL}/tcmspsearch.php",
        params={"qr": herb_en_name, "qsr": "herb_en_name", "token": token},
        timeout=30,
    )
    return _extract_grid_data(r.text)


def _extract_grid_data(html: str) -> list:
    """Extract the first Kendo Grid data array from page HTML."""
    # Find first occurrence of data: [{
    start = html.find('data: [{"')
    if start == -1:
        start = html.find('data:[{"')
    if start == -1:
        return []

    bracket_pos = html.index("[", start)
    decoder = json.JSONDecoder()
    try:
        data, _ = decoder.raw_decode(html[bracket_pos:])
        return data if isinstance(data, list) else []
    except Exception:
        return []


def _filter_compounds(raw: list) -> list:
    """Keep only compounds with OB >= threshold and DL >= threshold."""
    results = []
    for c in raw:
        try:
            ob = float(c.get("ob") or 0)
            dl = float(c.get("dl") or 0)
            if ob >= OB_THRESHOLD and dl >= DL_THRESHOLD:
                results.append({
                    "mol_name": (c.get("molecule_name") or "").strip(),
                    "OB": round(ob, 2),
                    "DL": round(dl, 2),
                    "MW": round(float(c.get("mw") or 0), 2),
                })
        except Exception:
            continue
    return results


def fetch_herb(session, token: str, herb_cn: str) -> list:
    """Fetch and filter compounds for one Chinese herb name."""
    # Step 1: search for the herb to get its English name
    candidates = _search_herbs(session, token, herb_cn)

    # Pick the best match (exact Chinese name first)
    herb_info = None
    for h in candidates:
        if h.get("herb_cn_name") == herb_cn:
            herb_info = h
            break
    if not herb_info and candidates:
        herb_info = candidates[0]

    # If Chinese search failed, try fallback English name
    if not herb_info and herb_cn in CN_TO_EN_FALLBACK:
        en_fallback = CN_TO_EN_FALLBACK[herb_cn]
        candidates2 = _search_herbs(session, token, en_fallback)
        if candidates2:
            herb_info = candidates2[0]

    if not herb_info:
        return []

    # Step 2: get compounds from detail page
    en_name = herb_info.get("herb_en_name", "")
    if not en_name:
        return []

    raw = _get_compounds(session, token, en_name)
    return _filter_compounds(raw)


def _save(db: dict):
    with open(BUILTIN_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, ensure_ascii=False, indent=2)


def main():
    # Load existing data
    if BUILTIN_PATH.exists():
        with open(BUILTIN_PATH, encoding="utf-8") as f:
            db = json.load(f)
        print(f"✅ 已有内置数据：{len(db)} 味中药")
    else:
        db = {}

    to_fetch = [h for h in HERBS if h not in db]
    print(f"📋 需要抓取：{len(to_fetch)} 味\n")

    if not to_fetch:
        print("全部已有数据，无需抓取。")
        return

    session = _make_session()
    token = _get_token(session)
    print(f"🔑 Token: {token}\n")

    for i, herb in enumerate(to_fetch):
        print(f"  [{i+1}/{len(to_fetch)}] {herb} ...", end=" ", flush=True)
        try:
            rows = fetch_herb(session, token, herb)
            if rows:
                db[herb] = rows
                print(f"✅ {len(rows)} 个活性成分")
            else:
                print("⚠️  无数据")
        except Exception as e:
            print(f"❌ {e}")

        # Save every 5 herbs
        if (i + 1) % 5 == 0:
            _save(db)
            print(f"    💾 已保存进度 ({len(db)} 味)")

        time.sleep(1.0)  # polite delay

    _save(db)
    print(f"\n🎉 完成！共 {len(db)} 味中药")
    print("\n下一步推送：")
    print("  git add data/tcmsp_builtin.json && git commit -m '补充TCMSP内置数据' && git push")


if __name__ == "__main__":
    main()
