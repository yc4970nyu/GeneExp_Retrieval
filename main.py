import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Title of the App
st.title("Gene Expression and Correlation Analysis Tool")
st.subheader("Visualize gene expression, explore correlations, and discover insights")

# Step 1: Load Dataset
@st.cache_data
def load_data():
    file_path = 'Anticancer_Drug_Treatment_DATA.txt'  # Ensure this file is in the same directory
    df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
    df.index.name = 'Gene'
    df.columns.name = 'Treatment'
    return df

df = load_data()

# Sidebar for input
st.sidebar.header("Gene Analysis")
gene_name = st.sidebar.text_input("Enter a gene name for individual analysis:")
st.sidebar.subheader("Select Up to 5 Genes for Group Analysis")
genes = []
for i in range(5):
    gene = st.sidebar.text_input(f"Gene {i + 1} (Optional):", key=f"gene_{i}")
    if gene:
        genes.append(gene)

# Main content area
if gene_name:
    if gene_name in df.index:
        # Gene Expression Visualization (Single Gene)
        st.subheader(f"📊 Expression Levels of {gene_name}")
        gene_data = df.loc[gene_name]
        
        # Line plot for single gene
        plt.figure(figsize=(12, 6))
        plt.plot(df.columns, gene_data, marker='o', label=f'Response of {gene_name}', color='blue')
        plt.title(f"Expression Levels of {gene_name} Across Treatments", fontsize=16)
        plt.xlabel("Treatments", fontsize=14)
        plt.ylabel("Expression Level", fontsize=14)
        plt.xticks(rotation=45, fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        st.pyplot(plt)

        # Top 10 Correlated Genes
        st.subheader(f"🔬 Top 10 Genes Correlated with {gene_name}")
        correlations = df.corrwith(df.loc[gene_name], axis=1)
        top_correlated_genes = correlations.drop(gene_name).sort_values(ascending=False).head(10)

        # Display the top 10 correlated genes as a table
        st.table(top_correlated_genes)

        # Bar plot of top 10 correlated genes
        plt.figure(figsize=(10, 6))
        sns.barplot(
            x=top_correlated_genes.values,
            y=top_correlated_genes.index,
            palette='coolwarm'
        )
        plt.title(f"Top 10 Genes Correlated with {gene_name}", fontsize=16)
        plt.xlabel("Pearson Correlation", fontsize=14)
        plt.ylabel("Genes", fontsize=14)
        plt.tight_layout()
        st.pyplot(plt)
    else:
        st.error(f"Gene '{gene_name}' not found in the dataset. Please try a different gene.")

if genes:
    valid_genes = [gene for gene in genes if gene in df.index]
    invalid_genes = [gene for gene in genes if gene not in df.index]

    # Feedback for invalid genes
    if invalid_genes:
        st.error(f"The following genes were not found: {', '.join(invalid_genes)}")

    if valid_genes:
        # Gene Expression Visualization (Multiple Genes)
        st.subheader("📊 Expression Levels of Selected Genes")
        plt.figure(figsize=(12, 8))
        palette = sns.color_palette("Set2", len(valid_genes))
        
        for idx, gene in enumerate(valid_genes):
            plt.plot(df.columns, df.loc[gene], marker='o', label=f'Response of {gene}', color=palette[idx])
        
        plt.title("Cardiotoxicity Response of Selected Genes to Treatments", fontsize=16)
        plt.xlabel("Treatments", fontsize=14)
        plt.ylabel("Cardiotoxicity Response", fontsize=14)
        plt.xticks(rotation=45, fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(title="Genes", fontsize=12, title_fontsize=14, loc='upper right')
        plt.grid(True)
        plt.tight_layout()
        st.pyplot(plt)

        # Correlation Matrix
        if len(valid_genes) >= 2:
            st.subheader("🔬 Correlation Matrix of Selected Genes")
            correlation_matrix = df.loc[valid_genes].T.corr()

            # Heatmap of correlation matrix
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
            plt.xticks(rotation=45, fontsize=12)
            plt.yticks(fontsize=12)
            plt.tight_layout()
            st.pyplot(plt)
else:
    if not gene_name:
        st.info("Please enter a gene name in the sidebar for individual analysis.")
    if not genes:
        st.info("Please select up to 5 genes in the sidebar for group analysis.")

# Optional: Display Raw Data
if st.checkbox("📂 Show raw dataset"):
    st.dataframe(df)


# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")


