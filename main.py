import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

# Section: Individual Gene Analysis
st.header("📊 Individual Gene High-Correlation Retrieval")
if gene_name:
    if gene_name in df.index:
        col1, col2 = st.columns([1, 1])
        
        # Gene Expression Line Plot
        with col1:
            st.subheader(f"Expression Levels of {gene_name}")
            plt.figure(figsize=(10, 6))
            plt.plot(df.columns, df.loc[gene_name], marker='o', label=f'{gene_name}', color='blue')
            plt.title(f"Expression Levels of {gene_name} Across Treatments", fontsize=16)
            plt.xlabel("Treatments", fontsize=14)
            plt.ylabel("Expression Level", fontsize=14)
            plt.xticks(rotation=45, fontsize=12)
            plt.legend(fontsize=12)
            plt.grid(True)
            plt.tight_layout()
            st.pyplot(plt)

        # Top 10 Correlated Genes
        with col2:
            st.subheader(f"Top 10 Correlated Genes with {gene_name}")
            correlations = df.corrwith(df.loc[gene_name], axis=1)
            top_correlated_genes = correlations.drop(gene_name).sort_values(ascending=False).head(10)
            st.table(top_correlated_genes)

            plt.figure(figsize=(8, 6))
            sns.barplot(
                x=top_correlated_genes.values,
                y=top_correlated_genes.index,
                palette="coolwarm"
            )
            plt.title(f"Top 10 Correlated Genes with {gene_name}", fontsize=16)
            plt.xlabel("Pearson Correlation", fontsize=14)
            plt.ylabel("Genes", fontsize=14)
            plt.tight_layout()
            st.pyplot(plt)
    else:
        st.error(f"Gene '{gene_name}' not found in the dataset.")

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

    # Validate yerr shape to avoid errors
    try:
        yerr = summary_stats["Std Expression"].values
        if len(yerr) != len(summary_stats["TT Code"]):
            raise ValueError("Invalid shape for yerr.")
    except Exception as e:
        st.warning(f"Error with error bars for {treatment}: {e}. Error bars will be skipped.")
        yerr = None

    # Bar plot with error bars
    plt.figure(figsize=(8, 6))
    sns.barplot(
        data=summary_stats,
        x="TT Code",
        y="Mean Expression",
        yerr=yerr,  # Validated error bars
        palette="gray",
        ci=None
    )

    # Overlay individual data points
    sns.stripplot(
        data=group_data,
        x="TT Code",
        y="Expression",
        color="black",
        size=6,
        jitter=True
    )

    # Add plot formatting
    plt.title(f"Expression Levels for {treatment}", fontsize=16)
    plt.xlabel("TT Code", fontsize=14)
    plt.ylabel("Expression Level", fontsize=14)
    plt.xticks(rotation=45, fontsize=12)
    plt.tight_layout()
    st.pyplot(plt)

# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")


