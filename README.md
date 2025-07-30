# LogP and TPSA Correlation Analysis

## 🔍 Overview
This project analyzes the correlation between **molecular lipophilicity (LogP)** and **Topological Polar Surface Area (TPSA)** using the `Lipophilicity.csv` dataset.

> 🧑‍🔬 This project was conducted independently by a first-year chemical engineering student  
> It demonstrates the use of cheminformatics and machine learning to explore and validate a chemical hypothesis through a computational analysis.

## 🎯 Objective
To validate the chemical hypothesis that:
> Molecules with higher TPSA values tend to have lower LogP values (i.e., more polar → less liposoluble), indicating a negative correlation between the two variables.

## 🙋 Feedback Wanted!
I'm actively improving this project and would love your feedback on:

- Code structure and readability  
- Analysis and visualization  
- Feature engineering ideas  
- Model evaluation and interpretation  

Feel free to open an [Issue](https://github.com/Kyeongbin-Tom/logp-tpsa-analysis/issues) or start a [Discussion](https://github.com/Kyeongbin-Tom/logp-tpsa-analysis/discussions) to share your thoughts!

## 🧭 Motivation
As a first-year chemical engineering student, I became curious about how TPSA (a measure of polarity) influences other molecular properties.  
I hypothesized that lipophilicity (LogP) would be a key factor, as it reflects a molecule’s solubility balance between water (polar) and oil (nonpolar) phases.  
This project represents my first attempt to apply data analysis and machine learning techniques to verify a chemical hypothesis.

## 🧠 Background
Lipophilicity (LogP) is a measure of how well a compound dissolves in fat relative to water. Topological Polar Surface Area (TPSA) quantifies the polar regions in the molecule, which are capable of forming hydrogen bonds. TPSA reflects polarity due to its relation to hydrogen bonding capacity. 

In general, highly polar molecules are more hydrophilic and less likely to dissolve in lipid environments.

This leads to the hypothesis that as TPSA increases, LogP tends to decrease - indicating a negative relationship.

## 🧪 Methods
- Data processing using pandas
- Visualization with matplotlib
- Correlation analysis using Python
- Machine learning: SVR (Support Vector Regression) to predict the relationship between LogP and TPSA

## 📁 Data Source
- **Lipophilicity.csv**: From the CAIP GitHub repository (CC0 1.0)  
  🔗 [https://github.com/doas1min/CAIP/blob/main/data/Lipophilicity.csv](https://github.com/doas1min/CAIP/blob/main/data/Lipophilicity.csv)
- **TPSA Reference**:  
  Ertl, P., Rohde, B., & Selzer, P. (2000).  
  *Topological polar surface area: A useful descriptor in 2D‑QSAR*.  
  Journal of Medicinal Chemistry, 43(20), 3714‑3717.

## ⚙️ Setup and Execution

### Option 1: Run on Google Colab

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kyeongbin-Tom/logp-tpsa-analysis/blob/main/LogP_TPSA_Analysis.ipynb)

### Option 2: Run locally

```bash
git clone https://github.com/Kyeongbin-Tom/logp-tpsa-analysis.git
cd logp-tpsa-analysis
pip install -r requirements.txt
jupyter notebook notebooks/LogP_TPSA_Analysis.ipynb
```

⚠️ Some cells use Colab-specific commands (e.g., %pip install, !wget) and may require modification when running locally.

## 📊 Key Findings
- Pearson correlation coefficient between TPSA and LogP: **-0.14**
- Scatter plots shows a weak inverse trend.
- Supports the hypothesis that higher polarity (TPSA) reduces lipophilicity (LogP).

## 📈 Visualizations

### TPSA vs LogP
This scatter plot shows the relationship between TPSA (Topological Polar Surface Area) and LogP (lipophilicity), suggesting a weak negative correlation.

![TPSA vs LogP](https://raw.githubusercontent.com/Kyeongbin-Tom/logp-tpsa-analysis/main/images/tpsa_vs_logp.png)

### Parity Plot
This plot compares actual and predicted LogP values across three models:  
**Linear Regression**, **SVR**, and **Random Forest Regressor**.

![Parity Plot](https://raw.githubusercontent.com/Kyeongbin-Tom/logp-tpsa-analysis/main/images/parity_plot_all_models.png)

### SHAP Summary Plot
The SHAP summary plot illustrates the impact of each feature on the predicted LogP values, using the Random Forest model.

![SHAP Summary Plot](https://raw.githubusercontent.com/Kyeongbin-Tom/logp-tpsa-analysis/main/images/shap_summary_plot.png)

## 📂 Folder Structure

logp-tpsa-analysis/
├── data/
│   └── Lipophilicity.csv              # Input dataset (LogP and SMILES)
├── images/
│   ├── parity_plot_all_models.png     # Parity plot: Actual vs Predicted LogP (Linear Regression, SVR, Random Forest)
│   ├── shap_summary_plot.png          # SHAP summary plot: Feature impact on LogP prediction
│   └── tpsa_vs_logp.png               # TPSA vs LogP scatter plot
├── notebooks/
│   └── LogP_TPSA_Analysis.ipynb       # Main Jupyter notebook for correlation and modeling
├── src/
│   ├── descriptors.py                 # RDKit-based molecular descriptor calculator
│   └── visualization.py               # Plotting functions for LogP and TPSA
├── requirements.txt                   # Required Python packages
└── README.md                          # Project overview and instructions

## 🧰 Python Modules

### 🔬 'src/descriptors.py'
- 'smiles_to_mol(smiles)': Converts a SMILES string into an RDKit Mol object.
- 'smiles_to_descriptors(smiles)': Calculates molecular descriptors including TPSA, molecular weight, number of hydrogen bond acceptors/donors, and number of rotatable bonds.  
  -> Utilizes RDKit's 'Descriptors' module.

### 📊 'src/visualization.py'
- plot_tpsa_vs_logp(df): Plots the scatter plot of between TPSA vs LogP.
- plot_distribution(df, col): Plots the distribution of specified variables using histplot and KDE curve.
  -> Uses 'matplotlib' and 'seaborn'.


## 🤔 Conclusion
The hypothesis was supported: TPSA appears to negatively correlate with LogP. The correlation was lower than initially expected, but machine learning models still helped demonstrate the relationship.

## 💡 Future Work
- Try deeplearning models for higher accuracy.
- Identify additional molecular descriptors that correlate strongly with LogP.
- Utilize SHAP or permutation for better better interpretation.
- Explore other algorithms beyond SVR, such as XGBoost
- Learn and apply Organic Chemistry concepts to deepen the analysis of relationship between LogP and various molecular descriptors.

## 📌 License

The dataset (Lipophilicity.csv) is distributed under the
CC0 1.0 Public Domain Dedication.
This project is for educational and non-commercial purposes only.


## 📝 Citation

Ertl, P., Rohde, B., & Selzer, P. (2000).
Topological polar surface area: A useful descriptor in 2D-QSAR.
Journal of Medicinal Chemistry, 43(20), 3714–3717.
