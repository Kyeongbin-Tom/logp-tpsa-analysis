# LogP and TPSA Correlation Analysis

## 🔍 Overview
이 프로젝트는 Lipophilicity 데이터셋을 사용해 **분자 지용성(LogP)** 과 **토폴로지컬 극성 표면적(TPSA)** 간의 상관관계를 분석합니다.
This project analyzes correlation between **molecules' fat solubility (LogP)** and **Topological Polar Surface Area (TPSA)** using Lipophilicity.csv dataset.

**hypothesis:**  
"더 높은 TPSA를 가진 분자는 일반적으로 더 낮은 LogP 값을 가지며 (즉, 더 극성 → 덜 지용성), 두 변수 사이에 반비례 상관관계가 있다."
"Molecules with higher TPSA values tend to have lower LogP values (i.e., more polar → less liposoluble), indicating an inverse relationship between the two variables."

## 📁 Data source
- Lipophilicity.csv: CAIP GitHub repository (CC0 1.0)  
  `https://github.com/doas1min/CAIP/blob/main/data/Lipophilicity.csv`
- TPSA reference: Ertl, P., Rohde, B., & Selzer, P. (2000). *Topological polar surface area: A useful descriptor in 2D‑QSAR*. J. Med. Chem. 43(20), 3714‑3717. :contentReference[oaicite:4]{index=4}

## 📊 Key Findings
- Pearson correlation coefficient between TPSA and LogP: **-0.14**
- Visualization shows a clear inverse relationship.
- Supports the hypothesis that higher polarity (TPSA) reduces lipophilicity (LogP).

- ## 🧪 How to Run
You can run the notebook via:
- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/doas1min/logp-tpsa-analysis/blob/main/notebooks/LogP_TPSA_Analysis.ipynb)
- Or clone this repo and run locally:
```bash
pip install pandas seaborn matplotlib
jupyter notebook notebooks/LogP_TPSA_Analysis.ipynb
