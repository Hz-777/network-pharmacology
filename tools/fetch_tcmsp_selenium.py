"""
用 Selenium 自动抓取 TCMSP 数据（本地运行，支持 JS 渲染）
用法：python3 tools/fetch_tcmsp_selenium.py
"""

import json, time, sys
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from webdriver_manager.chrome import ChromeDriverManager

BUILTIN_PATH = Path(__file__).parent.parent / "data" / "tcmsp_builtin.json"
TCMSP_URL    = "https://www.tcmsp-e.com/tcmsp.php"
OB_THRESHOLD = 30.0
DL_THRESHOLD = 0.18

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


def init_driver(headless: bool = True) -> webdriver.Chrome:
    options = webdriver.ChromeOptions()
    if headless:
        options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-blink-features=AutomationControlled")
    options.add_argument("--window-size=1400,900")
    options.add_experimental_option("excludeSwitches", ["enable-automation"])
    options.add_experimental_option("useAutomationExtension", False)

    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=options)
    driver.execute_script(
        "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
    )
    return driver


def _dismiss_popup(driver):
    """Remove the daily popup and any overlays blocking the page."""
    driver.execute_script("""
        // Set localStorage to suppress future popups
        try { localStorage.setItem('tcmsp_daily_popup_date', new Date().toISOString().slice(0,10)); } catch(e){}

        // Remove popup element
        var dp = document.getElementById('dp');
        if (dp) dp.remove();

        // Remove any modal backdrops
        document.querySelectorAll('.modal-backdrop, .b[aria-hidden]').forEach(function(el){ el.remove(); });

        // Restore body scroll
        document.body.style.overflow = '';
        document.body.style.paddingRight = '';
        document.documentElement.style.overflow = '';
    """)
    time.sleep(0.3)


def _js(driver, script, *args):
    return driver.execute_script(script, *args)


def search_herb(driver, herb: str) -> list:
    """搜索中药，提取活性成分列表。"""
    driver.get(TCMSP_URL)
    time.sleep(3)
    _dismiss_popup(driver)

    # Use JS to fill & submit the search form
    submitted = _js(driver, """
        // Find the search input
        var input = document.querySelector('input[name="qr"]') ||
                    document.querySelector('#qrInput') ||
                    document.querySelector('input[type="text"]');
        if (!input) return false;

        // Set search type to Chinese herb name
        var sel = document.querySelector('select[name="qsr"]');
        if (sel) {
            for (var i=0; i<sel.options.length; i++) {
                if (sel.options[i].value === 'herb_cn_name') {
                    sel.selectedIndex = i;
                    break;
                }
            }
        }

        // Set value and trigger events
        input.value = arguments[0];
        input.dispatchEvent(new Event('input', {bubbles:true}));
        input.dispatchEvent(new Event('change', {bubbles:true}));

        // Submit the form
        var form = input.closest('form') || document.querySelector('form');
        if (form) {
            form.submit();
            return true;
        }

        // Fallback: press Enter
        input.dispatchEvent(new KeyboardEvent('keydown', {key:'Enter', keyCode:13, bubbles:true}));
        input.dispatchEvent(new KeyboardEvent('keypress', {key:'Enter', keyCode:13, bubbles:true}));
        input.dispatchEvent(new KeyboardEvent('keyup', {key:'Enter', keyCode:13, bubbles:true}));
        return true;
    """, herb)

    if not submitted:
        return []

    time.sleep(4)
    _dismiss_popup(driver)

    # Look for herb result links and click the first matching one
    current_url = driver.current_url
    try:
        links = driver.find_elements(By.XPATH, "//table//a | //td//a")
        for link in links:
            try:
                text = link.text.strip()
                href = link.get_attribute("href") or ""
                if herb in text or ("tcmsp.php" in href and "herb_cn_name" not in href):
                    if link.is_displayed():
                        _js(driver, "arguments[0].click();", link)
                        time.sleep(4)
                        _dismiss_popup(driver)
                        break
            except Exception:
                continue
    except Exception:
        pass

    # Try clicking the "Ingredients" tab
    for tab_text in ["Ingredients", "Molecules", "成分", "活性成分"]:
        try:
            tabs = driver.find_elements(By.XPATH, f"//*[normalize-space(text())='{tab_text}']")
            for tab in tabs:
                if tab.is_displayed():
                    _js(driver, "arguments[0].click();", tab)
                    time.sleep(2)
                    break
        except Exception:
            continue

    return _extract_all_pages(driver)


def _extract_all_pages(driver) -> list:
    """Extract table data across all pages."""
    all_rows = []
    for _ in range(20):  # max 20 pages
        all_rows += _get_table_rows(driver)

        # Try to go to next page
        went_next = False
        for xpath in [
            "//a[contains(@class,'k-next-page')]",   # Kendo pager
            "//span[contains(@class,'k-next')]/..",
            "//a[@title='Go to the next page']",
            "//li[contains(@class,'next') and not(contains(@class,'disabled'))]/a",
            "//a[normalize-space(text())='Next']",
            "//a[normalize-space(text())='下一页']",
        ]:
            try:
                btns = driver.find_elements(By.XPATH, xpath)
                for btn in btns:
                    cls = btn.get_attribute("class") or ""
                    if btn.is_displayed() and "disabled" not in cls and "k-state-disabled" not in cls:
                        _js(driver, "arguments[0].click();", btn)
                        time.sleep(2)
                        went_next = True
                        break
                if went_next:
                    break
            except Exception:
                continue

        if not went_next:
            break

    # Filter by OB/DL thresholds
    results = []
    for row in all_rows:
        try:
            ob = float(str(row.get("OB", "0")).replace("%", "") or 0)
            dl = float(str(row.get("DL", "0")) or 0)
            if ob >= OB_THRESHOLD and dl >= DL_THRESHOLD:
                results.append({
                    "mol_name": row.get("mol_name", "").strip(),
                    "OB":  round(ob, 2),
                    "DL":  round(dl, 2),
                    "MW":  round(float(str(row.get("MW", "0")) or 0), 2),
                    "SMILES": row.get("SMILES", "").strip(),
                })
        except Exception:
            continue
    return results


def _get_table_rows(driver) -> list:
    """Extract rows from the current page's data table."""
    rows = []
    try:
        tables = driver.find_elements(By.TAG_NAME, "table")
        for table in tables:
            headers_els = table.find_elements(By.TAG_NAME, "th")
            headers = [h.text.strip() for h in headers_els]
            header_str = " ".join(headers).upper()

            if not any(k in header_str for k in ["OB", "DL", "MOLECULE", "MOL", "NAME"]):
                continue

            # Build column index
            col = {}
            for i, h in enumerate(headers):
                lo = h.lower().strip()
                if any(k in lo for k in ["mol_name", "molecule name", "compound", "name", "成分"]):
                    col.setdefault("mol_name", i)
                elif lo in ("ob", "ob (%)") or lo.startswith("ob"):
                    col["OB"] = i
                elif lo in ("dl",) or lo.startswith("dl"):
                    col["DL"] = i
                elif lo in ("mw", "mol_weight", "molecular weight"):
                    col["MW"] = i
                elif "smiles" in lo:
                    col["SMILES"] = i

            if "mol_name" not in col:
                col["mol_name"] = 0

            for tr in table.find_elements(By.TAG_NAME, "tr")[1:]:
                tds = tr.find_elements(By.TAG_NAME, "td")
                if not tds:
                    continue
                row = {}
                for key, idx in col.items():
                    if idx < len(tds):
                        row[key] = tds[idx].text.strip()
                if row:
                    rows.append(row)

            if rows:
                break  # found the right table

    except Exception:
        pass
    return rows


def _save(db: dict):
    with open(BUILTIN_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, ensure_ascii=False, indent=2)


def main():
    if BUILTIN_PATH.exists():
        with open(BUILTIN_PATH, encoding="utf-8") as f:
            db = json.load(f)
        print(f"✅ 已有内置数据：{len(db)} 味中药")
    else:
        db = {}

    to_fetch = [h for h in HERBS if h not in db]
    print(f"📋 需要抓取：{len(to_fetch)} 味")

    if not to_fetch:
        print("全部已有数据，无需抓取。")
        return

    print("🌐 启动浏览器（headless 模式）...\n")
    driver = init_driver(headless=True)

    try:
        for i, herb in enumerate(to_fetch):
            print(f"  [{i+1}/{len(to_fetch)}] {herb} ...", end=" ", flush=True)
            try:
                rows = search_herb(driver, herb)
                if rows:
                    db[herb] = rows
                    print(f"✅ {len(rows)} 个活性成分")
                else:
                    print("⚠️  无数据")
            except Exception as e:
                print(f"❌ 错误: {e}")

            if (i + 1) % 5 == 0:
                _save(db)
                print(f"    💾 已保存进度 ({len(db)} 味)")

            time.sleep(2)

    finally:
        driver.quit()
        _save(db)

    print(f"\n🎉 完成！共 {len(db)} 味中药")
    print("下一步推送：")
    print("  git add data/tcmsp_builtin.json && git commit -m '补充TCMSP内置数据' && git push")


if __name__ == "__main__":
    main()
