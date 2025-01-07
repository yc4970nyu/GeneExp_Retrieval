import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set page layout
st.set_page_config(layout="wide", page_title="Gene Expression and Treatment Analysis")

# Title
st.title("🔬 Gene Expression and Treatment Analysis Tool")
st.markdown("Analyze gene expression across treatments, explore correlations, and visualize grouped treatment responses with advanced tools.")

# Step 1: Load Dataset
@st.cache_data
def load_data():
    file_path = 'Anticancer_Drug_Treatment_DATA.txt'  # Ensure this file is in the same directory
    df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
    df.index.name = 'Gene'
    df.columns.name = 'Treatment'
    return df

df = load_data()

# Mapping TT codes to treatment types
treatment_mapping = {
    "TT-235": "DMSO", "TT-236": "DMSO", "TT-237": "DMSO", "TT-238": "DMSO", "TT-239": "DMSO", "TT-240": "DMSO",
    "TT-241": "Doxorubicin", "TT-242": "Doxorubicin",
    "TT-243": "Bortezomib", "TT-244": "Bortezomib",
    "TT-245": "Tegretol", "TT-246": "Tegretol",
    "TT-247": "WZ8040", "TT-248": "WZ8040",
    "TT-249": "PHA-767491", "TT-250": "PHA-767491",
    "TT-251": "BI 2536", "TT-252": "BI 2536",
    "TT-253": "Daunorubicin", "TT-254": "Daunorubicin",
    "TT-255": "Nitrendipine", "TT-256": "Nitrendipine",
    "TT-257": "Solifenacin", "TT-258": "Solifenacin"
}

# Sidebar for Inputs
st.sidebar.header("Input Options")
gene_name = st.sidebar.text_input("🔎 Enter a gene name for individual analysis:")
st.sidebar.subheader("Select Up to 5 Genes for Group Analysis")
genes = [st.sidebar.text_input(f"Gene {i + 1} (Optional):", key=f"gene_{i}") for i in range(5)]
genes = [gene for gene in genes if gene]

# Section: Treatment-Type Visualizations
st.header("🧪 Treatment-Type Visualizations")

# Group by treatment types
df_long = df.melt(var_name="TT Code", value_name="Expression")
df_long["Treatment"] = df_long["TT Code"].map(treatment_mapping)
treatment_groups = df_long.groupby("Treatment")

for treatment, group_data in treatment_groups:
    st.subheader(f"Treatment: {treatment}")
    
    # Calculate mean and standard deviation
    summary_stats = group_data.groupby("TT Code").agg({"Expression": ["mean", "std"]}).reset_index()
    summary_stats.columns = ["TT Code", "Mean Expression", "Std Expression"]

    # Ensure Std Expression has no NaN values
    summary_stats["Std Expression"] = summary_stats["Std Expression"].fillna(0)

    # Matplotlib bar plot with error bars
    plt.figure(figsize=(8, 6))
    x = np.arange(len(summary_stats["TT Code"]))  # X-axis positions
    means = summary_stats["Mean Expression"]
    errors = summary_stats["Std Expression"]

    plt.bar(x, means, yerr=errors, capsize=5, color='gray', alpha=0.7, edgecolor='black', label="Mean ± Std")
    plt.xticks(x, summary_stats["TT Code"], rotation=45, fontsize=12)
    plt.ylabel("Expression Level", fontsize=14)
    plt.xlabel("TT Code", fontsize=14)
    plt.title(f"Expression Levels for {treatment}", fontsize=16)
    plt.tight_layout()

    # Overlay individual data points
    for idx, tt_code in enumerate(summary_stats["TT Code"]):
        tt_data = group_data[group_data["TT Code"] == tt_code]["Expression"]
        jittered_x = np.random.normal(x[idx], 0.05, size=len(tt_data))  # Add jitter for visualization
        plt.scatter(jittered_x, tt_data, color='black', s=20, label="Data Points" if idx == 0 else "")

    plt.legend()
    st.pyplot(plt)

# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")


