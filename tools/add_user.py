#!/usr/bin/env python3
"""
添加新用户到 users.yaml
用法：python3 tools/add_user.py
"""

import sys
import os
import getpass
from pathlib import Path

# 确保从项目根目录运行
ROOT = Path(__file__).parent.parent
USERS_FILE = ROOT / "users.yaml"

try:
    import yaml
    import bcrypt
except ImportError:
    print("缺少依赖，请先运行：pip install PyYAML bcrypt")
    sys.exit(1)


def load_config() -> dict:
    if not USERS_FILE.exists():
        print(f"❌ 未找到 {USERS_FILE}，请先确认文件存在")
        sys.exit(1)
    with open(USERS_FILE, encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_config(config: dict):
    with open(USERS_FILE, "w", encoding="utf-8") as f:
        yaml.dump(config, f, allow_unicode=True, default_flow_style=False)


def list_users(config: dict):
    users = config.get("credentials", {}).get("usernames", {})
    if not users:
        print("（暂无用户）")
        return
    print(f"\n{'用户名':<15} {'显示名称':<12} {'邮箱':<30} {'角色'}")
    print("-" * 70)
    for uname, info in users.items():
        print(f"{uname:<15} {info.get('name',''):<12} {info.get('email',''):<30} {info.get('role','user')}")


def add_user(config: dict):
    print("\n── 添加新用户 ──")
    username = input("用户名（英文，用于登录）: ").strip()
    if not username:
        print("用户名不能为空")
        return

    users = config.setdefault("credentials", {}).setdefault("usernames", {})
    if username in users:
        print(f"⚠️  用户名 '{username}' 已存在")
        return

    name  = input("显示名称（中文可以）: ").strip() or username
    email = input("邮箱地址: ").strip()
    role  = input("角色 [user/admin]（默认 user）: ").strip() or "user"

    while True:
        pwd = getpass.getpass("登录密码: ")
        pwd2 = getpass.getpass("再次确认密码: ")
        if pwd == pwd2:
            break
        print("两次密码不一致，请重试")

    hashed = bcrypt.hashpw(pwd.encode(), bcrypt.gensalt()).decode()

    users[username] = {
        "name": name,
        "email": email,
        "password": hashed,
        "role": role,
    }
    save_config(config)
    print(f"\n✅ 用户 '{username}'（{name}）添加成功")


def delete_user(config: dict):
    users = config.get("credentials", {}).get("usernames", {})
    username = input("\n要删除的用户名: ").strip()
    if username not in users:
        print(f"⚠️  用户 '{username}' 不存在")
        return
    confirm = input(f"确认删除 '{username}'？(y/N): ").strip().lower()
    if confirm != "y":
        print("已取消")
        return
    del users[username]
    save_config(config)
    print(f"✅ 用户 '{username}' 已删除")


def change_password(config: dict):
    users = config.get("credentials", {}).get("usernames", {})
    username = input("\n要修改密码的用户名: ").strip()
    if username not in users:
        print(f"⚠️  用户 '{username}' 不存在")
        return
    while True:
        pwd = getpass.getpass("新密码: ")
        pwd2 = getpass.getpass("再次确认: ")
        if pwd == pwd2:
            break
        print("两次密码不一致，请重试")
    hashed = bcrypt.hashpw(pwd.encode(), bcrypt.gensalt()).decode()
    users[username]["password"] = hashed
    save_config(config)
    print(f"✅ 用户 '{username}' 密码已更新")


def main():
    print("=" * 40)
    print("  网络药理学平台 — 用户管理工具")
    print("=" * 40)

    config = load_config()

    while True:
        print("\n请选择操作：")
        print("  1. 查看所有用户")
        print("  2. 添加新用户")
        print("  3. 删除用户")
        print("  4. 修改密码")
        print("  0. 退出")
        choice = input("\n输入数字: ").strip()

        if choice == "1":
            list_users(config)
        elif choice == "2":
            add_user(config)
            config = load_config()  # reload
        elif choice == "3":
            list_users(config)
            delete_user(config)
            config = load_config()
        elif choice == "4":
            list_users(config)
            change_password(config)
            config = load_config()
        elif choice == "0":
            print("退出")
            break
        else:
            print("无效输入，请重试")


if __name__ == "__main__":
    main()
