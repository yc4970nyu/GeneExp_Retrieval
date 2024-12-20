import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Title of the App
st.title("Gene Expression and Correlation Analysis Tool")
st.subheader("Visualize gene expression and discover correlations")

# Step 1: Load Dataset
@st.cache
def load_data():
    file_path = 'Anticancer_Drug_Treatment_DATA.txt'  # Ensure this file is in the same directory
    df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
    df.index.name = 'Gene'
    df.columns.name = 'Treatment'
    return df

df = load_data()

# Sidebar for input
st.sidebar.header("Gene Analysis")
gene_name = st.sidebar.text_input("Enter a gene name for analysis:")

# Main content area
if gene_name:
    if gene_name in df.index:
        # Gene Expression Visualization
        st.subheader(f"📊 Expression Levels of {gene_name}")
        gene_data = df.loc[gene_name]
        
        # Line plot of gene expression
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

        # Correlation Analysis: Top 10 correlated genes
        st.subheader(f"🔬 Top 10 Genes Correlated with {gene_name}")
        correlations = df.corrwith(df.loc[gene_name], axis=1)
        top_correlated_genes = correlations.drop(gene_name).sort_values(ascending=False).head(10)

        # Display the top 10 correlated genes as a table
        st.table(top_correlated_genes)

        # Visualize the correlations in a bar plot
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
else:
    st.info("Please enter a gene name in the sidebar to start the analysis.")

# Optional: Display Raw Data
if st.checkbox("📂 Show raw dataset"):
    st.dataframe(df)


# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")


