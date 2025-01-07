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

# Section: New Feature - Improved Visualization with Box Plot
st.header("📊 Gene Expression Across All Treatment Types (Box Plot)")
if gene_name:
    if gene_name in df.index:
        # Filter data for the selected gene
        gene_data = df.loc[gene_name]

        # Normalize the gene data
        normalized_gene_data = gene_data / gene_data.max()

        # Prepare data for plotting
        gene_long = normalized_gene_data.reset_index()
        gene_long.columns = ["TT Code", "Normalized Expression"]
        gene_long["Treatment"] = gene_long["TT Code"].map(treatment_mapping)

        # Box plot with individual data points
        plt.figure(figsize=(12, 8))
        sns.boxplot(
            data=gene_long,
            x="Treatment",
            y="Normalized Expression",
            palette="Set2",
            showmeans=True,
            meanprops={"marker": "o", "markerfacecolor": "red", "markeredgecolor": "black"}
        )
        sns.stripplot(
            data=gene_long,
            x="Treatment",
            y="Normalized Expression",
            color="black",
            size=5,
            jitter=True,
            alpha=0.7
        )

        # Plot formatting
        plt.xlabel("Treatment Types", fontsize=14)
        plt.ylabel("Normalized Expression", fontsize=14)
        plt.title(f"Normalized Expression of {gene_name} Across All Treatment Types", fontsize=16)
        plt.xticks(rotation=45, fontsize=12)
        plt.tight_layout()
        st.pyplot(plt)
    else:
        st.error(f"Gene '{gene_name}' not found in the dataset.")
        
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
            plt.barh(top_correlated_genes.index, top_correlated_genes.values, color='skyblue')
            plt.title(f"Top 10 Correlated Genes with {gene_name}", fontsize=16)
            plt.xlabel("Pearson Correlation", fontsize=14)
            plt.ylabel("Genes", fontsize=14)
            plt.tight_layout()
            st.pyplot(plt)

# Section: Existing Features (Multi-Gene Analysis, Correlation Matrix)
st.header("📈 Group Gene Analysis")
if genes:
    valid_genes = [gene for gene in genes if gene in df.index]
    invalid_genes = [gene for gene in genes if gene not in df.index]

    if invalid_genes:
        st.error(f"The following genes were not found: {', '.join(invalid_genes)}")

    if valid_genes:
        # Multi-Gene Line Plot
        st.subheader("Expression Levels of Selected Genes")
        plt.figure(figsize=(12, 8))
        palette = plt.cm.Set2(np.linspace(0, 1, len(valid_genes)))
        for idx, gene in enumerate(valid_genes):
            plt.plot(df.columns, df.loc[gene], marker='o', label=f'{gene}', color=palette[idx])
        plt.title("Expression Levels of Selected Genes Across Treatments", fontsize=16)
        plt.xlabel("Treatments", fontsize=14)
        plt.ylabel("Expression Levels", fontsize=14)
        plt.xticks(rotation=45, fontsize=12)
        plt.legend(title="Genes", fontsize=12, title_fontsize=14)
        plt.grid(True)
        plt.tight_layout()
        st.pyplot(plt)

        # Correlation Matrix
        if len(valid_genes) > 1:
            st.subheader("Correlation Matrix of Selected Genes")
            correlation_matrix = df.loc[valid_genes].T.corr()
            plt.figure(figsize=(8, 6))
            plt.imshow(correlation_matrix, cmap="coolwarm", interpolation='nearest')
            plt.colorbar()
            plt.title("Correlation Matrix", fontsize=16)
            plt.xticks(range(len(valid_genes)), valid_genes, rotation=45, fontsize=10)
            plt.yticks(range(len(valid_genes)), valid_genes, fontsize=10)
            plt.tight_layout()
            st.pyplot(plt)

# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")




