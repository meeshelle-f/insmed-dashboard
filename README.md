# Protein Analysis Dashboard

This Streamlit web application provides interactive tools for exploring protein data, divided into two main sections: **Peptide Group Analysis** and **Deep Mutational Scanning (DMS)** Analysis.

## Features

### 1. Peptide Group Analysis

This section helps you understand relationships and similarities between different groups of peptides and how individual peptides relate to those groups.

* **Data Loading:** The application automatically loads peptide group data from a GitHub URL (`peptide_groups.csv`). It includes robust data validation and cleaning steps to ensure the presence of 'group' and 'peptide' columns with valid string entries. A dummy dataset is used if loading from the URL fails.

* **Reference Peptide Similarity:**
    * Input a **custom peptide sequence** (must be 9 characters long and use the allowed amino acid alphabet).
    * Calculates and visualizes the **similarity score distribution** between your reference peptide and all peptides within selected groups using a box plot. This helps you see how consistently similar a group's peptides are to your chosen reference.
    * Provides quantitative statistics (mean, median, standard deviation, IQR) to summarize the spread of similarity scores within each group.
* **Inter-Group Similarity:**
    * **Identifies the most similar pair of peptide groups** by calculating the average pairwise similarity between all unique combinations of groups.
    * Visualizes the **average similarity scores between all group pairs** in a bar chart, allowing you to easily compare group relationships.
    * Presents a **detailed pairwise similarity heatmap** for a selected group pair, showing the similarity score for every single peptide from one group against every peptide from another. This provides a granular view of inter-group relationships.

---

### 2. Deep Mutational Scanning (DMS) Analysis & Prediction Comparison

This section focuses on understanding the functional impact of single amino acid changes in a protein, displaying raw DMS data, and comparing it to computational predictions.

* **DMS Data Loading:** The application loads experimental DMS data (`dms.csv`), along with predicted scores from the ESM-IF1 (`if1.csv`) and ProteinMPNN (`mpnn.csv`) models from a GitHub URL. The data is preprocessed to extract original amino acid, position, and new amino acid from the 'mutant' column.
* **DMS Score Heatmap:** 
    * Generates a **heatmap of average DMS scores** across different protein positions and new amino acids. This helps identify positions and mutations that significantly affect protein fitness or function.
    * Allows you to **select a specific position range** to focus the heatmap.
    * Provides **dynamic sizing and conditional annotation** for optimal visualization.
    * Highlights the **easiest and hardest positions to mutate** and the most/least tolerated new amino acids based on average DMS scores.

* **DMS vs. Prediction Comparison:**
    * **Calculates and displays Pearson correlation coefficients** between experimental DMS scores and predicted scores from the IF1 and MPNN models, giving an overall measure of model accuracy.
    * **Identifies which model shows higher similarity** to the experimental DMS data.
    * Offers a **detailed comparison** for a selected model (IF1 or MPNN):
        * A **scatter plot** visualizes the relationship between experimental and predicted scores, with a diagonal line indicating perfect agreement.
        * Tables show the **top 10 mutations with the highest agreement** (lowest absolute difference between experimental and predicted scores) and the **top 10 mutations with the lowest agreement** (highest absolute difference), helping pinpoint where models perform well or poorly.

## Data Sources

The application fetches data from the following GitHub URLs:

* **Peptide Group Data:** `https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/peptides/peptide_groups.csv`
* **Deep Mutational Scanning (DMS) Data:** `https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/dms.csv`
* **ESM-IF1 Prediction Data:** `https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/if1.csv`
* **ProtMPNN Prediction Data:** `https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/mpnn.csv`

The application includes robust error handling and will use dummy data for relevant sections if the online data sources are unavailable or corrupted.

## How to Run

To run this Streamlit application locally, follow these steps:

1.  **Clone the repository** 
2.  **Install the required libraries:**
    ```bash
    pip install streamlit pandas numpy matplotlib seaborn
    ```
3.  **Navigate to the directory** containing your Streamlit app file (e.g., `app.py`) in your terminal.
4.  **Run the application:**
    ```bash
    streamlit run app.py
    ```
    This will open the application in your default web browser.

## Technologies Used

* **Streamlit:** For creating interactive web applications with Python.
* **Pandas:** For data manipulation and analysis.
* **NumPy:** For numerical operations.
* **Matplotlib & Seaborn:** For data visualization.
* **`itertools` & `re`:** For efficient iteration and regular expression operations.

---