# LogP and TPSA Correlation Analysis

## 🔍 Overview
This project analyzes the correlation between **molecules' fat solubility (LogP)** and **Topological Polar Surface Area (TPSA)** using the `Lipophilicity.csv` dataset.

**Hypothesis:**  
"Molecules with higher TPSA values tend to have lower LogP values (i.e., more polar → less liposoluble), indicating an inverse relationship between the two variables."

---

## 📁 Data Source
- **Lipophilicity.csv**: From the CAIP GitHub repository (CC0 1.0)  
  🔗 [https://github.com/doas1min/CAIP/blob/main/data/Lipophilicity.csv](https://github.com/doas1min/CAIP/blob/main/data/Lipophilicity.csv)
- **TPSA Reference**:  
  Ertl, P., Rohde, B., & Selzer, P. (2000).  
  *Topological polar surface area: A useful descriptor in 2D‑QSAR*.  
  Journal of Medicinal Chemistry, 43(20), 3714‑3717.

---

## 📊 Key Findings
- Pearson correlation coefficient between TPSA and LogP: **-0.14**
- Visualization shows a clear inverse relationship.
- Supports the hypothesis that higher polarity (TPSA) reduces lipophilicity (LogP).

---

## 🧪 How to Run

You can run the notebook via:

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kyeongbin-Tom/logp-tpsa-analysis/blob/main/LogP_TPSA_Analysis.ipynb)

- Or clone this repo and run locally:
  ```bash
  pip install pandas seaborn matplotlib
  jupyter notebook notebooks/LogP_TPSA_Analysis.ipynb

⚠️ This notebook is designed for Google Colab.
Some commands (e.g., %pip install, !wget) may not work in local Jupyter environments without modification.

## 📂 Folder Structure
logp-tpsa-analysis/
├── Lipophilicity.csv # Input dataset
├── LogP_TPSA_Analysis.ipynb # Main analysis notebook
└── README.md # Project overview


## 📌 License

The dataset (Lipophilicity.csv) is distributed under the
CC0 1.0 Public Domain Dedication.
This project is for educational and non-commercial purposes only.


## 📝 Citation

Ertl, P., Rohde, B., & Selzer, P. (2000).
Topological polar surface area: A useful descriptor in 2D-QSAR.
Journal of Medicinal Chemistry, 43(20), 3714–3717.
