#!/bin/bash
# 启动网络药理学分析平台
cd "$(dirname "$0")"
export PATH="$PATH:/Users/mac/Library/Python/3.9/bin"
python3 -m streamlit run app.py --server.port 8501 --server.headless false
