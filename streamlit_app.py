import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import re 

# --- Peptide Similarity Calculation ---
def calculate_peptide_similarity(peptide_1: str, peptide_2: str) -> float:
    """
    Calculates the proportion of identical amino acids between two peptides. 
    Function assumes input peptides are of equal length. 

    Args: 
        peptide_1 (str): The first peptide sequence string. 
        peptide_2 (str): The second peptide sequence string. 

    Returns:
        float: A similarity score between 0 and 1, representing the proportion of matching amino acids at identical positons. 
        Returns 0 if peptides are of difference lengths or are empty 
    """
    # 1. Type Coercion, Input Cleaning 
    peptide_1, peptide_2 = str(peptide_1).strip(), str(peptide_2).strip()

    # 2. Input Validation: handle empty or differing length strings
    if not peptide_1 or not peptide_2:
        return 0.0

    if len(peptide_1) != len(peptide_2):
        return 0.0

    # 3. Calculate 
    num_matched = 0
    for i in range(len(peptide_1)):
        if peptide_1[i] == peptide_2[i]:
            num_matched += 1
            
    return num_matched / len(peptide_1)


# --- Calculate similarity scores for a reference peptide ---
@st.cache_data #streamlit decorator 
def calculate_similarity_scores(reference_peptide: str, df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates similarity scores between a reference peptide and peptides in a DataFrame. 

    Scores are calculated using calculate_peptide_similarity. If a peptide in the
    DataFrame is invalid (empty or non-string leading to error), its score will be NaN.

    Args:
        reference_peptide (str): The peptide sequence to compare against. 
        df (pd.DataFrame): DataFrame containing a 'peptide' column. 

    Returns:
        pd.DataFrame(): A new DataFrame with an added 'similarity_score' column. 
                        Scores are NaN if the reference peptide is invalid or comparison fails. 
    """

    df_copy = df.copy()
    df_copy['similarity_score'] = float('nan')
    
    # 1. Clean and validate reference peptide once 
    reference_peptide = str(reference_peptide).strip()
    if not reference_peptide:
        return df_copy

    for index, row in df_copy.iterrows():
        peptide_x = str(row['peptide']).strip()
        if not peptide_x:
            continue
        try:
            score = calculate_peptide_similarity(reference_peptide, peptide_x)
            df_copy.at[index, 'similarity_score'] = score
        except (ValueError, TypeError):
            continue
    return df_copy


# --- Determine most similar peptide groups ---
@st.cache_data # decorator caches the function's output
def determine_most_similar_groups(df: pd.DataFrame) -> dict:
    """
    Determines the most similar pair of peptide groups within a DataFrame.

    This function calculates the average pairwise similarity score betwen all possible unique combinations of groups in the DataFrame's 'group' column.
    For each pair of groups, similarity is computed between every peptide from the first group and every peptide from the second group.
    The average of these individual similarities represents the overall similarity between the group pair. 
    The group with the highest overall average similarity score is returned. 

    Args:
        df (pd.DataFrame): DataFrame containing a 'group' and 'peptide' column. 
                        -'group': Identifies the group each peptide belongs to. 
                        -'peptide': Peptide sequence string. 
    Returns:
        group_pair (tuple): A tuple `(group1_name, group2_name)` representing the pair of groups with the highest avg similarity score. 
        average_similarity_score (float): The similairty score. Returns -1.0 if `group_pair` is `None`. 
        all_pair_scores (dict): Keys are `(group1_name, group2_name)` and values are their calculated avg similarity scores. 

    Notes:
    - `calculate_peptide_similarity()` calculates similarity between 2 peptides

    Example: 
    # df = pd.DataFrame({'group': ['G1', 'G1', 'G2', 'G2'],'peptide': ['AB', 'AB', 'AA', 'BB']})
    # calculate_peptide_similarity(AB,AA) = .5,  calculate_peptide_similarity(AB,BB) = .5
    #result = {'group_pair': ('G1', 'G2'), 'average_similarity_score': 0.5, 'all_pair_scores': {('G1', 'G2'): 0.5}}

    """

    unique_groups = df['group'].dropna().unique().tolist()

    if len(unique_groups) < 2: #early exit if there are fewer than two unique groups 
        return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}

    group_pairs = list(itertools.combinations(unique_groups, 2))

    all_group_pair_scores = {}

    for group1_name, group2_name in group_pairs:
        peptides_group1 = df[df['group'] == group1_name]['peptide'].dropna().tolist()
        peptides_group2 = df[df['group'] == group2_name]['peptide'].dropna().tolist()
        pairwise_scores_for_current_pair = [] #store all individual similarity scores 

        if not peptides_group1 or not peptides_group2: #early exit if either groups in pair is invalid
            all_group_pair_scores[(group1_name, group2_name)] = 0.0
            continue

        for p1 in peptides_group1:
            if not isinstance(p1, str) or not p1.strip(): continue
            for p2 in peptides_group2:
                if not isinstance(p2, str) or not p2.strip(): continue
                try:
                    score = calculate_peptide_similarity(p1.strip(), p2.strip())
                    pairwise_scores_for_current_pair.append(score)
                except (ValueError, TypeError):
                    continue
        
        if pairwise_scores_for_current_pair:
            avg_score = sum(pairwise_scores_for_current_pair) / len(pairwise_scores_for_current_pair)
            all_group_pair_scores[(group1_name, group2_name)] = avg_score
        else:
            all_group_pair_scores[(group1_name, group2_name)] = 0.0

    if not all_group_pair_scores:
        return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}
        
    most_similar_pair = max(all_group_pair_scores.items(), key=lambda item: item[1]) #pull max item from dictionary 

    return {
        'group_pair': most_similar_pair[0],
        'average_similarity_score': most_similar_pair[1],
        'all_pair_scores': all_group_pair_scores
    }


# --- Data Loading Functions ---
@st.cache_data
def load_peptide_groups_data():
    """
    Loads peptide group data from GitHub URL. 
    Performs essential data validation and cleaning steps, ensuring 'group' and 'peptide' columns are present, are of string type and contain non-empty values. 
    In case of loading errors or missing columns, dummy DataFrame is generated. 
    """
    try:
        df_loaded = pd.read_csv("https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/peptides/peptide_groups.csv")
        if 'group' not in df_loaded.columns or 'peptide' not in df_loaded.columns:
            st.error("Error: 'peptide_groups.csv' must contain 'group' and 'peptide' columns.")
            return pd.DataFrame()
        df_loaded['peptide'] = df_loaded['peptide'].astype(str).str.strip()
        df_loaded = df_loaded[df_loaded['peptide'] != ''].copy()
        df_loaded['group'] = df_loaded['group'].astype(str).str.strip()
        df_loaded = df_loaded[df_loaded['group'] != ''].copy()
        return df_loaded
    except Exception as e:
        st.error(f"Error loading default 'peptide_groups.csv' from URL: {e}")
        st.warning("Using a dummy 'Peptide Groups' DataFrame for demonstration due to loading failure.")
        data = {
            'group': ['group_A']*50 + ['group_B']*50 + ['group_C']*50,
            'peptide': [
                'AYLDVFNKF', 'VYIDVFNWD', 'AALDVFNKF', 'PYLDVFNKW', 'AYLDFNKKF',
            ]*10 + [
                'VYIDVFNWD', 'GYIDVFNWD', 'TYLDVFNWS', 'VYIDVENWA', 'AYLDVFNKE',
            ]*10 + [
                'PPLDVFNKA', 'QYLDVFNKC', 'RYIDVFNKD', 'AYLDVFNKG', 'AYLDVFNKH'
            ]*10
        }
        import random
        import string
        def generate_random_peptide(length=9):
            return ''.join(random.choice(string.ascii_uppercase) for _ in range(length))
        for i in range(len(data['peptide'])):
            if random.random() < 0.2:
                data['peptide'][i] = generate_random_peptide()
        return pd.DataFrame(data)

# Helper for parsing mutant string
mutant_pattern = re.compile(r'([A-Z])(\d+)([A-Z])') #[A-Z] captures single uppercase letter, (\d+) captures one or more digits 

def parse_mutant_column(df, mutant_col='mutant'):
    
    """
    Parses a specified mutant column in a DataFrame into separate 'original_aa', 'position', and 'new_aa' columns.
    """
    # Convert column to string type first, handling potential NaNs
    df[mutant_col] = df[mutant_col].astype(str)
    extracted_data = df[mutant_col].str.extract(mutant_pattern, expand=True)
    extracted_data.columns = ['original_aa', 'position_str', 'new_aa']

    #Convert 'position_str' to numeric
    extracted_data['position'] = pd.to_numeric(extracted_data['position_str'], errors='coerce').astype('Int64')
    extracted_data = extracted_data.drop(columns=['position_str'])
    # Concatenate the new columns back to the original DataFrame
    df = pd.concat([df, extracted_data], axis=1)
    return df 



#Loading 3 proteingym datasets, highthroughput measures functional effects of single mutations in a protein. 

@st.cache_data
def load_dms_data():
    """
    Loads and preprocesses Deep Mutational Scanning (DMS) experimental data.
    """
    try:
        dms_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/dms.csv')
        if 'mutant' not in dms_df.columns or 'DMS_score' not in dms_df.columns:
            st.error("Error: 'dms.csv' must contain 'mutant' and 'DMS_score' columns.")
            return pd.DataFrame()
        
        dms_df = parse_mutant_column(dms_df.copy())
        dms_df.dropna(subset=['original_aa', 'position', 'new_aa', 'DMS_score'], inplace=True)
        dms_df['position'] = dms_df['position'].astype(int)
        
        return dms_df
    except Exception as e:
        st.error(f"Error loading 'dms.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
        return pd.DataFrame()

@st.cache_data
def load_if1_data():
    try:
        if1_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/if1.csv')
        if 'mutant' not in if1_df.columns or 'esmif1_ll' not in if1_df.columns:
            st.error("Error: 'if1.csv' must contain 'mutant' and 'esmif1_ll' columns.")
            return pd.DataFrame()
        
        if1_df = parse_mutant_column(if1_df.copy())
        if1_df.dropna(subset=['original_aa', 'position', 'new_aa', 'esmif1_ll'], inplace=True)
        if1_df['position'] = if1_df['position'].astype(int)
        
        return if1_df
    except Exception as e:
        st.error(f"Error loading 'if1.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
        return pd.DataFrame()

@st.cache_data
def load_mpnn_data():
    try:
        mpnn_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/mpnn.csv')
        if 'mutant' not in mpnn_df.columns or 'pmpnn_ll' not in mpnn_df.columns:
            st.error("Error: 'mpnn.csv' must contain 'mutant' and 'pmpnn_ll' columns.")
            return pd.DataFrame()
        
        mpnn_df = parse_mutant_column(mpnn_df.copy())
        mpnn_df.dropna(subset=['original_aa', 'position', 'new_aa', 'pmpnn_ll'], inplace=True)
        mpnn_df['position'] = mpnn_df['position'].astype(int)
        
        return mpnn_df
    except Exception as e:
        st.error(f"Error loading 'mpnn.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
        return pd.DataFrame()




# --- Streamlit UI App Layout ---
# 1. Peptide Group Analysis: understanding relationships (similarities) between different groups of peptides and how individual peptides relate to those groups. 
# 2. Deep Mutational Scanning (DMS) Analysis: 

st.set_page_config(layout="wide") # allows broswer page to take up more horizontal space, ideal for plots and tables. 
st.title("Protein Analysis Dashboard")
st.write("Explore peptide group similarities, deep mutational scanning (DMS) data, and compare predictions.")

# --- Creates Tabbed Interface for 2 Main Types of Protein Analyses ---
tab_peptide, tab_dms = st.tabs(["Peptide Group Analysis", "Deep Mutational Scanning (DMS) Analysis"])

# --- Peptide Group Analysis Tab ---
with tab_peptide:
    st.header("Peptide Group Similarity Dashboard")
    st.write("This dashboard allows you to explore peptide group similarities based on sequence matching.")

    # --- Data Source Selection for Peptide Analysis ---
    data_source_option = st.radio(
        "Choose data source for Peptide Analysis:",
        ("Use existing 'Peptide Groups' data", ), # Tuple with single element
        key="peptide_data_source_radio"
    )

    df_peptide = None
    if data_source_option == "Use existing 'Peptide Groups' data":
        df_peptide = load_peptide_groups_data()
        if not df_peptide.empty:
            st.success("Default 'Peptide Groups' data loaded successfully!")
            
    # --- Main Peptide Dashboard Logic (Conditional on Data Loading) ---
    if df_peptide is not None and not df_peptide.empty:
        # Peptide Data Overview moved here to only appear with this tab
        st.sidebar.header("Peptide Data Overview")
        st.sidebar.write(f"**Total Peptides:** {len(df_peptide)}")
        st.sidebar.write(f"**Number of Groups:** {df_peptide['group'].nunique()}")
        
        st.subheader("Peptide Data Preview")
        st.dataframe(df_peptide.head())


        # --- Section 1: Analyze Reference Peptide Similarity ---
        st.header("1. Reference Peptide Similarity Analysis")
        st.write("Determine how similar a chosen peptide is to all peptides within each group.")

        st.sidebar.subheader("Reference Peptide Configuration")

        selected_reference_peptide = None 
        is_reference_peptide_valid = False # Flag for validation

        all_peptides_str = "".join(df_peptide['peptide'].dropna().tolist())
        peptide_alphabet = sorted(list(set(char.upper() for char in all_peptides_str if char.isalpha())))
        
        if peptide_alphabet:
            st.sidebar.info(f"**Allowed peptide characters:** {' '.join(peptide_alphabet)}")
        else:
            st.sidebar.warning("Could not determine a peptide alphabet from loaded peptide data. Please ensure 'peptide' column contains valid sequences.")

        # Default to custom peptide input
        # Note: max_peptide_len here is just for the input box, not a strict validation constraint
        max_peptide_len_from_data = df_peptide['peptide'].str.len().max() if not df_peptide.empty else 20
        # Set a reasonable max for user input if data is empty or provides no max len
        display_max_chars = max(9, max_peptide_len_from_data) 


        custom_peptide_input = st.sidebar.text_input(
            "Enter your custom peptide (must be 9 characters long and use allowed characters, e.g., 'AYLDVFNKF'):",
            value='AYLDVFNKF', # Default value for custom input
            max_chars=display_max_chars, # Allow user to type beyond 9 but validation will catch it
            key="custom_ref_peptide_input"
        )
        custom_peptide_input_clean = custom_peptide_input.strip().upper()

        if custom_peptide_input_clean:
            if not custom_peptide_input_clean.isalpha():
                st.sidebar.error("Custom peptide must contain only alphabetical characters (A-Z, a-z).")
                selected_reference_peptide = None
            elif len(custom_peptide_input_clean) != 9: # New length validation
                st.sidebar.error("Custom peptide must be exactly 9 characters long.")
                selected_reference_peptide = None
            else:
                is_valid_chars = True
                invalid_chars = []
                if peptide_alphabet: 
                    for char in custom_peptide_input_clean:
                        if char not in peptide_alphabet:
                            is_valid_chars = False
                            invalid_chars.append(char)
                
                if not is_valid_chars and peptide_alphabet: 
                    st.sidebar.error(f"Invalid characters detected: {' '.join(set(invalid_chars))}. Only these are allowed: {' '.join(peptide_alphabet)}")
                    selected_reference_peptide = None
                elif not is_valid_chars and not peptide_alphabet: 
                    st.sidebar.error("Invalid characters detected. Only alphabetical characters are allowed.")
                    selected_reference_peptide = None
                else:
                    selected_reference_peptide = custom_peptide_input_clean
                    is_reference_peptide_valid = True # Set flag to True only if all checks pass
        else:
            st.sidebar.info("Please enter a peptide sequence.")
            selected_reference_peptide = None

        #identify all unique groups from pre-loaded data
        all_groups = df_peptide['group'].dropna().unique().tolist()
        if not all_groups:
            st.warning("No peptide groups found in the loaded data to analyze.")
        elif selected_reference_peptide: # execute IF there are available peptide groups AND a selected_reference_peptide was determined
            #multiselect widget
            selected_groups_for_ref = st.multiselect(
                "Select Groups to Analyze (for reference peptide):",
                options=all_groups,
                default=all_groups,
                key="ref_peptide_group_multiselect"
            )

            if selected_groups_for_ref:
                # select appropriate portion of df matching group of interest 
                df_filtered_for_ref = df_peptide[df_peptide['group'].isin(selected_groups_for_ref)]
                
                # calculate similarity scores, drop invalid scores
                df_with_scores = calculate_similarity_scores(selected_reference_peptide, df_filtered_for_ref)
                valid_scores_df = df_with_scores.dropna(subset=['similarity_score'])

                if not valid_scores_df.empty:
                    st.subheader(f"Similarity Score Distribution to '{selected_reference_peptide}'")
                    fig1, ax1 = plt.subplots(figsize=(12, 7))
                    sns.boxplot(x='group', y='similarity_score', data=valid_scores_df, ax=ax1, palette='viridis')
                    ax1.set_title(f'Similarity Score Distribution to "{selected_reference_peptide}" Across Selected Groups')
                    ax1.set_xlabel('Peptide Group')
                    ax1.set_ylabel(f'Similarity Score to "{selected_reference_peptide}"')
                    ax1.set_ylim(0, 1)
                    ax1.grid(axis='y', linestyle='--', alpha=0.7) #axis = y sets horizontal grid lines 
                    st.pyplot(fig1) # takes matplotlib object and creates it as image in st app 

                    st.subheader("Quantitative Spread within Each Group:")
                    #where x represents the Series of similarity_score for each individual group 
                    group_spread_stats = valid_scores_df.groupby('group')['similarity_score'].agg(
                        mean_score='mean',
                        median_score='median',
                        std_dev_score='std',
                        iqr_score=lambda x: x.quantile(0.75) - x.quantile(0.25),
                        count='count'
                    ).dropna()
                    st.dataframe(group_spread_stats.style.format("{:.4f}")) #display as interactive table in st app
                else:
                    st.info(f"No valid similarity scores found for '{selected_reference_peptide}' within the selected groups. This might be due to peptide length mismatch or group content.")
            else:
                st.info("Please select at least one group to analyze against the reference peptide.")
        else:
            st.info("Please enter a valid reference peptide (9 characters, allowed characters) to proceed with this analysis section.")


        # --- Section 2: Explore Inter-Group Similarity ---
        # This entire section only displays if the reference peptide is valid
        if is_reference_peptide_valid:
            st.header("2. Inter-Group Similarity Analysis")
            st.write("Identify and visualize how similar different peptide groups are to each other.")

            group_similarity_results = determine_most_similar_groups(df_peptide.copy())

            if group_similarity_results['group_pair'] and group_similarity_results['all_pair_scores']:
                sorted_pair_scores_list = sorted(group_similarity_results['all_pair_scores'].items(), key=lambda item: item[1], reverse=True) # decsending order
                # convert sorted list of group pairs into a DataFrame
                pair_scores_df = pd.DataFrame(
                    [(f"{pair[0]} vs {pair[1]}", score) for pair, score in sorted_pair_scores_list],
                    columns=['Group Pair', 'Average Similarity Score']
                )

                st.subheader("Average Similarity Scores Between All Group Pairs:")
                # Inter-Group Similarity, Barchart
                num_top_pairs = st.slider("Show Top N Most Similar Pairs:", 1, len(pair_scores_df), min(5, len(pair_scores_df)), key="num_top_pairs_slider")
                
                fig2, ax2 = plt.subplots(figsize=(12, max(6, num_top_pairs * 0.8)))
                sns.barplot(x='Average Similarity Score', y='Group Pair', data=pair_scores_df.head(num_top_pairs), palette='coolwarm', ax=ax2)
                ax2.set_title(f'Top {num_top_pairs} Most Similar Peptide Group Pairs')
                ax2.set_xlabel('Average Similarity Score (0-1)')
                ax2.set_ylabel('Peptide Group Pair')
                ax2.set_xlim(0, 1)
                ax2.grid(axis='x', linestyle='--', alpha=0.7)
                st.pyplot(fig2)

                st.success(f"**The single most similar pair is: {group_similarity_results['group_pair'][0]} and {group_similarity_results['group_pair'][1]}** (Average Score: {group_similarity_results['average_similarity_score']:.4f})")

                st.subheader("Detailed Pairwise Similarity Heatmap")
                
                # Unlike previous bar chart, this heatmap provides a detailed, peptide-level view of similarity between two specific groups. 
                # Also allows zoom-in of similarity score for every single peptide from one group against every single peptide from another group.
                
                available_pairs_for_heatmap = [f"{pair[0]} vs {pair[1]}" for pair, score in sorted_pair_scores_list]
                selected_heatmap_pair_str = st.selectbox(
                    "Select a Group Pair for Detailed Heatmap:",
                    options=available_pairs_for_heatmap,
                    index=0,
                    key="heatmap_pair_selectbox"
                )

                selected_group1_name, selected_group2_name = selected_heatmap_pair_str.split(' vs ')

                peptides_g1 = df_peptide[df_peptide['group'] == selected_group1_name]['peptide'].dropna().tolist()
                peptides_g2 = df_peptide[df_peptide['group'] == selected_group2_name]['peptide'].dropna().tolist()

                # allows user to determine number of peptides included in heatmap
                # default is 20 peptides for initial display, sets dynamic maximum that is capped at 50
                max_heatmap_peptides = st.slider(
                    "Limit # of peptides for heatmap (adjust for clarity):",
                    min_value=5, max_value=min(50, max(len(peptides_g1), len(peptides_g2), 10)),
                    value=min(20, max(len(peptides_g1), len(peptides_g2), 10)),
                    key="heatmap_peptide_limit"
                )
                peptides_g1_sampled = peptides_g1[:max_heatmap_peptides]
                peptides_g2_sampled = peptides_g2[:max_heatmap_peptides]

                # will hold the pairwise similarity scores
                #row and col names are the peptides from g1, g2  
                similarity_matrix_inter_group = pd.DataFrame(
                    np.nan,
                    index=peptides_g1_sampled,
                    columns=peptides_g2_sampled,
                    dtype=float
                )

                for p1 in peptides_g1_sampled:
                    for p2 in peptides_g2_sampled:
                        if not isinstance(p1, str) or not p1.strip() or not isinstance(p2, str) or not p2.strip():
                            similarity_matrix_inter_group.loc[p1, p2] = np.nan
                            continue
                        try:
                            score = calculate_peptide_similarity(p1, p2)
                            similarity_matrix_inter_group.loc[p1, p2] = score
                        except (ValueError, TypeError):
                            similarity_matrix_inter_group.loc[p1, p2] = np.nan
                            continue

                # Inter-Group Similarity, Heatmap

                # ensure the matrix isn't empty, and not entirely filled with NaN values 
                if not similarity_matrix_inter_group.empty and not similarity_matrix_inter_group.isna().all().all():
                    fig4, ax4 = plt.subplots(figsize=(max(10, len(peptides_g2_sampled)*0.8), max(8, len(peptides_g1_sampled)*0.8)))
                    sns.heatmap(
                        similarity_matrix_inter_group,
                        cmap='viridis',
                        annot=True if (len(peptides_g1_sampled) < 15 and len(peptides_g2_sampled) < 15) else False,
                        fmt=".2f",
                        cbar_kws={'label': 'Similarity Score (0-1)'},
                        ax=ax4
                    )
                    # annot=True if... Displays numerical similarity score inside each cell of the heatmap (annot=True) only if both groups have fewer than 15 peptides.  
                    # .2f formats annotations to 2 decimal places if annot is True. 
                   
                    ax4.set_title(f'Pairwise Similarity Between {selected_group1_name} and {selected_group2_name}')
                    ax4.set_xlabel(f'Peptides in {selected_group2_name}')
                    ax4.set_ylabel(f'Peptides in {selected_group1_name}')
                    plt.tight_layout()  # adjusts plot parameters, preventing labels from overlapping if tight layout 
                    st.pyplot(fig4)
                else:
                    st.info("No valid data or peptides for this heatmap. Adjust selections or check data content (e.g., peptide strings might be empty).")
            else:
                st.warning("Not enough distinct peptide groups to perform inter-group similarity analysis (need at least 2), or no valid peptide data loaded.")
        else:
            st.info("Please enter a valid reference peptide (9 characters long, with allowed characters) to enable 'Inter-Group Similarity Analysis'.")

    else:
        st.info("Please load peptide data using the options above to enable analysis.")




# --- Deep Mutational Scanning (DMS) Analysis & Prediction Comparison Tab ---
# focuses on understanding the functional impace of single amino acid changes in a protein. 
# It's split into two main sections: displaying the raw DMS data as a heatmap and then comparing it to computational predictions. 

# This DMS tab provides a comprehensive way to visualize the impact of mutations from experimental data
# and evaluate how well computational models can predict these experimental outcomes, helping identify which models are more reliable and for which types of mutations.


with tab_dms:
    # Data Loading and Initial Setup
    st.header("Deep Mutational Scanning (DMS) Analysis & Prediction Comparison")
    st.write("Visualize DMS scores as a heatmap to identify positions and mutations affecting protein fitness/function, and compare them with predicted scores.")

    dms_df_main = load_dms_data() # Main DMS data for this tab
    if1_df_main = load_if1_data() # Main IF1 data for this tab
    mpnn_df_main = load_mpnn_data() # Main MPNN data for this tab

    # --- Section 1: DMS Score Heatmap and Interpretation ---
    st.subheader("1. DMS Score Heatmap Analysis")
    if not dms_df_main.empty: #only run if DMS data was loaded successfully 
        st.success("DMS data loaded successfully!")
        st.subheader("DMS Data Preview")
        st.dataframe(dms_df_main.head())

        st.subheader("DMS Heatmap Configuration")

        # Heatmap Config and Data Prep
        all_positions = sorted(dms_df_main['position'].unique()) # extract all unique protein positions
        all_new_aas = sorted(dms_df_main['new_aa'].unique()) # extract all 'new' mutations 

        if not all_positions or not all_new_aas:
            st.warning("Not enough data to generate heatmap (missing positions or new amino acids after parsing).")
        else:
            min_pos, max_pos = min(all_positions), max(all_positions)
            selected_pos_range = st.slider(
                "Select Position Range for Heatmap:",
                min_value=min_pos,
                max_value=max_pos,
                value=(min_pos, max_pos),
                key="dms_pos_range_slider"
            )
 
            # filtered_dms_df: the DataFrame is filtered based on user's selected position range. 
            filtered_dms_df = dms_df_main[(dms_df_main['position'] >= selected_pos_range[0]) & 
                                     (dms_df_main['position'] <= selected_pos_range[1])]
            # Reshapes df into a matrix suitable for heatmap
            # each row represents a new_aa and each col represents position in the protein 
            # if there are multiple DMS scores for the same mutation at a specific position, this calculates their average 
            heatmap_data = filtered_dms_df.pivot_table(
                values='DMS_score',
                index='new_aa',
                columns='position',
                aggfunc='mean'
            ).fillna(np.nan)

            # Reindex to ensure all AAs and positions are present for consistent plotting
            heatmap_data = heatmap_data.reindex(index=all_new_aas, columns=range(min_pos, max_pos + 1), fill_value=np.nan)
            
            st.subheader("DMS Score Heatmap: Average Score per Mutation")
            st.write("Colors represent the average DMS score for mutations. Brighter (yellow) indicates higher scores, suggesting mutations are more tolerated at that position to that new amino acid.")

            #Plotting the heatmap 
            valid_scores_flat = heatmap_data.stack().dropna()
            if not valid_scores_flat.empty:
                default_vmin = valid_scores_flat.min()
                default_vmax = valid_scores_flat.max()
            else:
                default_vmin = 0.0
                default_vmax = 1.0

            #st.number_input widget allows manual control of heatmap color scale  
            col1, col2 = st.columns(2)
            with col1:
                heatmap_vmin = st.number_input("Heatmap Min Score (vmin):", value=default_vmin, key="heatmap_vmin", format="%.2f")
            with col2:
                heatmap_vmax = st.number_input("Heatmap Max Score (vmax):", value=default_vmax, key="heatmap_vmax", format="%.2f")

            if not heatmap_data.empty and not heatmap_data.isna().all().all():
                # --- Improved Heatmap Plotting ---
                num_aas = len(heatmap_data.index)
                num_positions = len(heatmap_data.columns)

                # Dynamically adjust figure size for better visibility. based on how many AA's and positions are being displayed.
                fig_dms, ax_dms = plt.subplots(figsize=(max(15, num_positions * 0.7), max(10, num_aas * 0.5))) 
                
                # Conditional annotation based on number of cells, only display numerical score within each haeatmap cell if there are relatively few cells (<15 and <12)
                annot_val = True if (num_aas < 15 and num_positions < 20) else False
                annot_kws = {"size": 8} if annot_val else {} 

                sns.heatmap(
                    heatmap_data,
                    cmap='viridis', # Good default, consider 'plasma', 'cividis', 'inferno' if needed
                    annot=annot_val, # Use conditional annotation
                    fmt=".2f",
                    cbar_kws={'label': 'Average DMS Score'},
                    vmin=heatmap_vmin,
                    vmax=heatmap_vmax,
                    ax=ax_dms,
                    annot_kws=annot_kws # Apply annotation font size
                )
                #labels and layout 
                ax_dms.set_title(f'Average DMS Score Heatmap (Positions {selected_pos_range[0]} - {selected_pos_range[1]})', fontsize=16)
                ax_dms.set_xlabel('Protein Position', fontsize=14)
                ax_dms.set_ylabel('New Amino Acid', fontsize=14)
                ax_dms.tick_params(axis='x', labelsize=10, rotation=90) # Rotate x-ticks for many positions
                ax_dms.tick_params(axis='y', labelsize=10, rotation=0) # Keep y-ticks horizontal
                plt.tight_layout() # Adjust layout to prevent labels from overlapping
                st.pyplot(fig_dms)
            else:
                st.info("No data available for the selected position range to generate the heatmap. Adjust the range or check data availability in the preview.")

            st.subheader("Interpretation: Easiest & Hardest to Mutate Positions")
            
            # --- Colorful Underlines ---
            # used to embed a small HTML snippet 
            st.markdown(
                """
                <p>
                (Assuming a higher DMS score indicates <span style="border-bottom: 3px solid yellow;">easier/more tolerated</span> mutation, 
                and a lower score indicates <span style="border-bottom: 3px solid purple;">harder/less tolerated</span> mutation)
                </p>
                """, 
                unsafe_allow_html=True
            )
            # avg_score_per_position: Calculates the avg DMS score for each protein position across all mutations at that position.
            avg_score_per_position = filtered_dms_df.groupby('position')['DMS_score'].mean().sort_values(ascending=False)
            # Calculates the average DMS score for each new amino acid across all positions. 
            avg_score_per_new_aa = filtered_dms_df.groupby('new_aa')['DMS_score'].mean().sort_values(ascending=False)

            # allows display of 'Positions by Overall Mutability (Average Score at Position)' & 'New Amino Acids by General Tolerability (Average Score to New AA)'
            col_pos, col_aa = st.columns(2)

            with col_pos:
                if not avg_score_per_position.empty:
                    st.write("#### Positions by Overall Mutability (Average Score at Position)")
                    st.write(f"**Easiest to Mutate Positions (Highest Average Score):**")
                    st.dataframe(avg_score_per_position.head(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
                    st.write(f"**Hardest to Mutate Positions (Lowest Average Score):**")
                    st.dataframe(avg_score_per_position.tail(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
                else:
                    st.info("No position-wise averages available for the selected range.")

            with col_aa:
                if not avg_score_per_new_aa.empty:
                    st.write("#### New Amino Acids by General Tolerability (Average Score to New AA)")
                    st.write(f"**Most Tolerated Mutations (to AA with highest avg score):**")
                    st.dataframe(avg_score_per_new_aa.head(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
                    st.write(f"**Least Tolerated Mutations (to AA with lowest avg score):**")
                    st.dataframe(avg_score_per_new_aa.tail(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
                else:
                    st.info("No new amino acid-wise averages available for the selected range.")
    else:
        st.info("DMS data is not loaded or is empty for heatmap analysis. Please check the URL or parsing in the code.")


    # --- Section 2: Prediction Comparison ---
    st.markdown("---") # Add a separator for visual clarity
    st.header("2. DMS vs. Prediction Comparison")
    st.write("Compare experimental Deep Mutational Scanning (DMS) scores with predicted scores from different models (esmif1_ll, pmpnn_ll).")

    if not dms_df_main.empty and not if1_df_main.empty and not mpnn_df_main.empty:
        st.success("All DMS, IF1, and MPNN prediction data loaded successfully for comparison!")

        # Merge for IF1 comparison
        comparison_df_if1 = pd.merge(dms_df_main[['mutant', 'DMS_score']], 
                                     if1_df_main[['mutant', 'esmif1_ll']], 
                                     on='mutant', 
                                     how='inner')
        
        # Merge for MPNN comparison
        comparison_df_mpnn = pd.merge(dms_df_main[['mutant', 'DMS_score']], 
                                      mpnn_df_main[['mutant', 'pmpnn_ll']], 
                                      on='mutant', 
                                      how='inner')
        # how = 'inner', only mutations present in both DataFrames will be included
        
        if comparison_df_if1.empty and comparison_df_mpnn.empty:
            st.warning("No common mutations found between DMS and either prediction dataset. Please check 'mutant' column formatting and data content.")
        else:
            st.subheader("Overall Model Similarity to DMS Data")

            #calculate correlation coefficients between experimental DMS_score and predicted scores (esmif1_ll or pmpnn_ll).

            correlation_if1 = 0.0
            if not comparison_df_if1.empty:
                correlation_if1 = comparison_df_if1['DMS_score'].corr(comparison_df_if1['esmif1_ll'])

            correlation_mpnn = 0.0
            if not comparison_df_mpnn.empty:
                correlation_mpnn = comparison_df_mpnn['DMS_score'].corr(comparison_df_mpnn['pmpnn_ll'])

            col_corr1, col_corr2 = st.columns(2)
            with col_corr1:
                st.metric(label="Pearson Correlation (DMS vs. IF1)", value=f"{correlation_if1:.4f}")
                if comparison_df_if1.empty:
                    st.info("No common mutations for IF1 correlation.")
            with col_corr2:
                # st.metric is fancy widget great for dispalying KPIs, making them stand out 
                st.metric(label="Pearson Correlation (DMS vs. MPNN)", value=f"{correlation_mpnn:.4f}")
                if comparison_df_mpnn.empty:
                    st.info("No common mutations for MPNN correlation.")
            
            st.markdown("---") # Visual separator

            if not comparison_df_if1.empty or not comparison_df_mpnn.empty:
                if correlation_if1 > correlation_mpnn:
                    st.success(f"**The IF1 model (Pearson R: {correlation_if1:.4f}) is more similar to DMS data than the MPNN model (Pearson R: {correlation_mpnn:.4f}).**")
                elif correlation_mpnn > correlation_if1:
                    st.success(f"**The MPNN model (Pearson R: {correlation_mpnn:.4f}) is more similar to DMS data than the IF1 model (Pearson R: {correlation_if1:.4f}).**")
                else:
                    st.info("The IF1 and MPNN models have very similar correlation values with DMS data.")
            else:
                st.warning("Not enough data to determine which model is more similar (no common mutations with DMS).")

            # Detailed Comparison (Scatter Plot and Tables)
            st.markdown("---") # Another visual separator
            st.subheader("Detailed Comparison for Selected Model")

            # allow user to choose which model to view 
            model_to_view = st.selectbox(
                "Select Prediction Model to View Details:",
                options=['IF1 (esmif1_ll)', 'MPNN (pmpnn_ll)'],
                key="model_comparison_selector"
            )

            if model_to_view == 'IF1 (esmif1_ll)':
                current_comparison_df = comparison_df_if1
                predicted_score_col = 'esmif1_ll'
                model_name = 'IF1'
            else: # MPNN (pmpnn_ll)
                current_comparison_df = comparison_df_mpnn
                predicted_score_col = 'pmpnn_ll'
                model_name = 'MPNN'

            if not current_comparison_df.empty:
                # Scatter Plot for selected model
                st.write(f"#### Scatter Plot: Experimental vs. Predicted Scores for {model_name}")
                st.write("Each point represents a mutation. Closeness to the diagonal line (y=x) indicates better agreement.")

                fig_scatter, ax_scatter = plt.subplots(figsize=(10, 8))
                # generate scatter plot with the DMS_score (experimental) on x-axis and predicted_score_col (model prediction) on the y-axis
                sns.scatterplot(x='DMS_score', y=predicted_score_col, data=current_comparison_df, ax=ax_scatter, alpha=0.6)
                
                max_val = max(current_comparison_df['DMS_score'].max(), current_comparison_df[predicted_score_col].max())
                min_val = min(current_comparison_df['DMS_score'].min(), current_comparison_df[predicted_score_col].min())
                #draw dashed red line where y=x, points closer to this line indicate better agreement between experimental and predicted scores 
                ax_scatter.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', label='Perfect Agreement')

                ax_scatter.set_title(f'DMS Score vs. {predicted_score_col} ({model_name} Predicted Score)')
                ax_scatter.set_xlabel('DMS Score (Experimental)')
                ax_scatter.set_ylabel(f'{predicted_score_col} (Predicted)')
                ax_scatter.legend()
                ax_scatter.grid(True, linestyle='--', alpha=0.7)
                st.pyplot(fig_scatter)

                # calculate absolute difference between experimental and predicted scores for each mutation 
                current_comparison_df['absolute_difference'] = abs(current_comparison_df['DMS_score'] - current_comparison_df[predicted_score_col])

                # Highest Agreement (Lowest Difference) for selected model
                st.write(f"#### Mutations with Highest Agreement ({model_name} Model)")
                st.write("These mutations have the smallest absolute difference between experimental and predicted scores.")
                highest_agreement_mutations = current_comparison_df.sort_values(by='absolute_difference', ascending=True).head(10)
                st.dataframe(highest_agreement_mutations[['mutant', 'DMS_score', predicted_score_col, 'absolute_difference']].style.format({"DMS_score": "{:.4f}", predicted_score_col: "{:.4f}", "absolute_difference": "{:.4f}"}))

                # Lowest Agreement (Highest Difference) for selected model
                st.write(f"#### Mutations with Lowest Agreement ({model_name} Model)")
                st.write("These mutations have the largest absolute difference between experimental and predicted scores, indicating significant disagreement.")
                lowest_agreement_mutations = current_comparison_df.sort_values(by='absolute_difference', ascending=False).head(10)
                st.dataframe(lowest_agreement_mutations[['mutant', 'DMS_score', predicted_score_col, 'absolute_difference']].style.format({"DMS_score": "{:.4f}", predicted_score_col: "{:.4f}", "absolute_difference": "{:.4f}"}))
            else:
                st.info(f"No common mutations found for {model_name} to perform detailed comparison.")

    elif dms_df_main.empty:
        st.warning("DMS data not available for prediction comparison. Please ensure DMS data loads successfully.")
    elif if1_df_main.empty:
        st.warning("Prediction data (if1.csv) not loaded. Cannot perform prediction comparison.")
    elif mpnn_df_main.empty:
        st.warning("Prediction data (mpnn.csv) not loaded. Cannot perform prediction comparison.")
    else:
        st.info("Loading all data for prediction comparison...")