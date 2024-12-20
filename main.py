import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set page layout
st.set_page_config(layout="wide", page_title="Gene Expression and Correlation Analysis")

# Title
st.title("🔬 Gene Expression and Correlation Analysis Tool")
st.markdown("Analyze gene expression across treatments, explore correlations, and uncover insights with advanced visualizations.")

# Step 1: Load Dataset
@st.cache_data
def load_data():
    file_path = 'Anticancer_Drug_Treatment_DATA.txt'  # Ensure this file is in the same directory
    df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
    df.index.name = 'Gene'
    df.columns.name = 'Treatment'
    return df

df = load_data()

# Sidebar for Inputs
st.sidebar.header("Input Options")
gene_name = st.sidebar.text_input("🔎 Enter a gene name for individual analysis:")
st.sidebar.subheader("Select Up to 5 Genes for Group Analysis")
genes = [st.sidebar.text_input(f"Gene {i + 1} (Optional):", key=f"gene_{i}") for i in range(5)]
genes = [gene for gene in genes if gene]

# Section: Individual Gene Analysis
st.header("📊 Individual Gene High-Corr Retrieval")
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

# Section: Group Gene Analysis
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
        palette = sns.color_palette("Set2", len(valid_genes))
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
            sns.heatmap(
                correlation_matrix,
                annot=True,
                cmap="coolwarm",
                fmt=".2f",
                cbar=True,
                xticklabels=valid_genes,
                yticklabels=valid_genes,
                linewidths=0.5
            )
            plt.title("Correlation Matrix of Selected Genes", fontsize=16)
            plt.tight_layout()
            st.pyplot(plt)

# Section: Creative Visualizations
st.header("🎨 Point-wise Visualizations")

# Heatmap of Gene-Treatment Responses
st.subheader("Heatmap of Selected Gene Responses Across Treatments")
if genes:
    selected_gene_data = df.loc[valid_genes]
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        selected_gene_data,
        cmap="coolwarm",
        annot=False,
        cbar=True,
        xticklabels=True,
        yticklabels=True
    )
    plt.title("Heatmap of Selected Gene Responses Across Treatments", fontsize=16)
    plt.xlabel("Treatments", fontsize=14)
    plt.ylabel("Genes", fontsize=14)
    plt.tight_layout()
    st.pyplot(plt)

# Interactive Scatter Plot
st.subheader("Interactive Gene-Treatment Scatter Plot")
if gene_name and gene_name in df.index:
    gene_data = df.loc[gene_name]
    treatment = st.selectbox("Select a Treatment for Scatter Plot:", df.columns)
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x=gene_data.index,
        y=gene_data.values,
        color="blue",
        s=100
    )
    plt.axvline(x=treatment, color="red", linestyle="--", label=f"Selected: {treatment}")
    plt.title(f"Expression Levels of {gene_name} Across Treatments", fontsize=16)
    plt.xlabel("Treatments", fontsize=14)
    plt.ylabel("Expression Level", fontsize=14)
    plt.legend()
    plt.grid(True)
    st.pyplot(plt)



# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")


