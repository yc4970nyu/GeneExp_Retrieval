import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Title of the App
st.title("Gene Analysis Tool")
st.subheader("Visualize gene expression and correlation analysis")

# Step 1: Load Dataset
@st.cache
def load_data():
    file_path = 'Anticancer_Drug_Treatment_DATA.txt'  # Ensure this file is in the same directory
    df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
    df.index.name = 'Gene'
    df.columns.name = 'Treatment'
    return df

df = load_data()

# Step 2: User Input for Gene Selection
st.subheader("Select Up to 5 Genes")
genes = []
for i in range(5):
    gene = st.text_input(f"Gene {i + 1} (Leave empty if not using):", "")
    if gene:
        genes.append(gene)

# Step 3: Visualization of Selected Genes
if genes:
    valid_genes = [gene for gene in genes if gene in df.index]
    invalid_genes = [gene for gene in genes if gene not in df.index]

    if invalid_genes:
        st.error(f"The following genes were not found: {', '.join(invalid_genes)}")
    
    if valid_genes:
        st.subheader("Expression Levels of Selected Genes")
        
        # Plot expression levels for valid genes
        plt.figure(figsize=(12, 8))
        for gene in valid_genes:
            plt.plot(df.columns, df.loc[gene], marker='o', label=f'Response of {gene}')
        
        plt.title("Cardiotoxicity Response of Selected Genes to Treatments", fontsize=16)
        plt.xlabel("Treatments", fontsize=14)
        plt.ylabel("Cardiotoxicity Response", fontsize=14)
        plt.xticks(rotation=45, fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(title="Genes", fontsize=12, title_fontsize=14, loc='upper right')
        plt.grid(True)
        plt.tight_layout()
        st.pyplot(plt)

# Step 4: Correlation Analysis
if len(valid_genes) >= 2:
    st.subheader("Correlation Analysis")
    
    # Compute correlation matrix
    correlation_matrix = df.loc[valid_genes].T.corr()
    
    # Visualize the correlation matrix as a heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', cbar=True, xticklabels=valid_genes, yticklabels=valid_genes)
    plt.title("Correlation Matrix of Selected Genes", fontsize=16)
    plt.xticks(rotation=45, fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    st.pyplot(plt)

# Step 5: Display Raw Data (Optional)
if st.checkbox("Show raw data"):
    st.write(df)

