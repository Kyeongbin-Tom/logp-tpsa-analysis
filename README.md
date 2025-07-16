# LogP and TPSA Correlation Analysis

## 🔍 Overview
This project analyzes the correlation between **molecular lipophilicity (LogP)** and **Topological Polar Surface Area (TPSA)** using the `Lipophilicity.csv` dataset.

## 🎯 Objective
To validate the chemical hypothesis that:
> Molecules with higher TPSA values tend to have lower LogP values (i.e., more polar → less liposoluble), indicating a negative correlation between the two variables.

## 🧪 Methods
- Data processing using pandas
- Visualization with matplotlib
- Correlation analysis by python
- Machine learning: SVR (Support Vector Regression) to predict the relationship between LogP and TPSA

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
- Scatter plots reveal a slight inverse trend.
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

## 📈 Visualizations

### TPSA vs LogP
![TPSA vs LogP](images/tpsa_vs_logp.png)

### SVR Actual vs Predicted
![SVR vs Actual](images/svr_vs_actual.png)


## 📂 Folder Structure
logp-tpsa-analysis/
├── Lipophilicity.csv # Input dataset
├── LogP_TPSA_Analysis.ipynb # Main analysis notebook
└── README.md # Project overview

## 🤔 Conclusion
The hypothesis was supported. TPSA is a significant predictor for LogP. The correlation was lower than initially expected, Machine learning models effectively demonstrated the relationship.

## 💡 Future Work
- Try deeplearning models for higher accuracy.
- Find other elements which has more relative with LogP

## 📌 License

The dataset (Lipophilicity.csv) is distributed under the
CC0 1.0 Public Domain Dedication.
This project is for educational and non-commercial purposes only.


## 📝 Citation

Ertl, P., Rohde, B., & Selzer, P. (2000).
Topological polar surface area: A useful descriptor in 2D-QSAR.
Journal of Medicinal Chemistry, 43(20), 3714–3717.
