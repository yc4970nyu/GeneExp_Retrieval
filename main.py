import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set Streamlit page configuration
st.set_page_config(page_title="Gene Analysis Tool", layout="wide", initial_sidebar_state="expanded")

# Custom styles for the app
st.markdown(
    """
    <style>
    .reportview-container {
        background: linear-gradient(135deg, #ffffff, #e3f2fd);
    }
    .sidebar .sidebar-content {
        background: #f7f7f7;
    }
    h1, h2, h3 {
        color: #2c3e50;
        font-family: 'Arial Black', sans-serif;
    }
    .css-1d391kg {
        font-family: 'Verdana', sans-serif;
        font-size: 16px;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# Title and header
st.title("Gene Expression Analysis Tool")
st.markdown("""
Welcome to the **Gene Expression Visualization Tool** owned by Koba Lab!  
Use this application to explore gene expression levels across treatments and perform advanced correlation analysis with a touch of artistry.  
""")

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
st.sidebar.header("Select Genes for Analysis")
genes = []
for i in range(5):
    gene = st.sidebar.text_input(f"Gene {i + 1} (Optional):", "")
    if gene:
        genes.append(gene)

# Main content area
if genes:
    valid_genes = [gene for gene in genes if gene in df.index]
    invalid_genes = [gene for gene in genes if gene not in df.index]

    # Error feedback for invalid genes
    if invalid_genes:
        st.error(f"The following genes were not found: {', '.join(invalid_genes)}")
    
    if valid_genes:
        # Gene Expression Visualization
        st.subheader("📊 Expression Levels of Selected Genes")
        plt.figure(figsize=(12, 8))
        palette = sns.color_palette("Set2", len(valid_genes))
        
        for idx, gene in enumerate(valid_genes):
            plt.plot(df.columns, df.loc[gene], marker='o', label=f'Response of {gene}', color=palette[idx])
        
        plt.title("Cardiotoxicity Response of Selected Genes to Treatments", fontsize=18, color="#34495e", weight="bold")
        plt.xlabel("Treatments", fontsize=14, color="#2c3e50")
        plt.ylabel("Cardiotoxicity Response", fontsize=14, color="#2c3e50")
        plt.xticks(rotation=45, fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(title="Genes", fontsize=12, title_fontsize=14, loc='upper right')
        plt.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        st.pyplot(plt)

        # Correlation Analysis
        if len(valid_genes) >= 2:
            st.subheader("🔬 Correlation Analysis of Selected Genes")
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
            plt.title("Correlation Matrix of Selected Genes", fontsize=18, color="#34495e", weight="bold")
            plt.xticks(rotation=45, fontsize=12)
            plt.yticks(fontsize=12)
            plt.tight_layout()
            st.pyplot(plt)
else:
    st.info("🎯 Please enter at least one gene in the sidebar to start the analysis.")

# Display Raw Data (Optional)
with st.expander("📂 View Raw Dataset"):
    st.dataframe(df)

# Footer
st.markdown("""
---
Designed by Louis Cui and supervised by Dr. Satoru Kobayashi.  
""")


