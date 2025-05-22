# import streamlit as st
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import itertools
# import io

# # --- Peptide Similarity Calculation Function (re-use from previous code) ---
# def calculate_peptide_similarity(peptide_1: str, peptide_2: str) -> float:
#     # Ensure peptides are processed as strings and stripped
#     peptide_1 = str(peptide_1).strip()
#     peptide_2 = str(peptide_2).strip()

#     if not peptide_1 or not peptide_2: # Handle empty strings
#         return 0.0

#     if len(peptide_1) != len(peptide_2):
#         # Decide how to handle different lengths: return 0, raise error, or calculate partial?
#         # For simplicity in this dashboard, we return 0.0 for different lengths.
#         return 0.0
            
#     num_matched = 0
#     for i in range(len(peptide_1)):
#         if peptide_1[i] == peptide_2[i]:
#             num_matched += 1
            
#     return num_matched / len(peptide_1)

# # --- Function to calculate similarity scores for a reference peptide (re-use) ---
# @st.cache_data # Cache this function to speed up reruns
# def calculate_similarity_scores(reference_peptide: str, df: pd.DataFrame) -> pd.DataFrame:
#     df_copy = df.copy()
#     df_copy['similarity_score'] = float('nan')
    
#     # Ensure reference_peptide is clean
#     reference_peptide = str(reference_peptide).strip()
#     if not reference_peptide: # If reference peptide is empty, no scores can be calculated
#         return df_copy

#     for index, row in df_copy.iterrows():
#         peptide_x = str(row['peptide']).strip() # Ensure it's a string and strip whitespace
#         if not peptide_x: # Skip empty peptide strings in the dataframe
#             continue
#         try:
#             score = calculate_peptide_similarity(reference_peptide, peptide_x)
#             df_copy.at[index, 'similarity_score'] = score
#         except (ValueError, TypeError): # Catch errors from calculate_peptide_similarity
#             continue
#     return df_copy

# # --- Function to determine most similar groups (re-use) ---
# @st.cache_data # Cache this function as it can be computationally intensive
# def determine_most_similar_groups(df: pd.DataFrame) -> dict:
#     unique_groups = df['group'].dropna().unique().tolist()
#     if len(unique_groups) < 2:
#         return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}

#     group_pairs = list(itertools.combinations(unique_groups, 2))
#     all_group_pair_scores = {}

#     for group1_name, group2_name in group_pairs:
#         peptides_group1 = df[df['group'] == group1_name]['peptide'].dropna().tolist()
#         peptides_group2 = df[df['group'] == group2_name]['peptide'].dropna().tolist()
#         pairwise_scores_for_current_pair = []

#         if not peptides_group1 or not peptides_group2:
#             all_group_pair_scores[(group1_name, group2_name)] = 0.0
#             continue

#         for p1 in peptides_group1:
#             if not isinstance(p1, str) or not p1.strip(): continue # Skip non-strings or empty
#             for p2 in peptides_group2:
#                 if not isinstance(p2, str) or not p2.strip(): continue # Skip non-strings or empty
#                 try:
#                     score = calculate_peptide_similarity(p1.strip(), p2.strip())
#                     pairwise_scores_for_current_pair.append(score)
#                 except (ValueError, TypeError):
#                     continue
        
#         if pairwise_scores_for_current_pair:
#             avg_score = sum(pairwise_scores_for_current_pair) / len(pairwise_scores_for_current_pair)
#             all_group_pair_scores[(group1_name, group2_name)] = avg_score
#         else:
#             all_group_pair_scores[(group1_name, group2_name)] = 0.0

#     if not all_group_pair_scores:
#         return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}
        
#     most_similar_pair = max(all_group_pair_scores.items(), key=lambda item: item[1])

#     return {
#         'group_pair': most_similar_pair[0],
#         'average_similarity_score': most_similar_pair[1],
#         'all_pair_scores': all_group_pair_scores
#     }

# # --- Data Loading (from your previous code, now with caching) ---
# @st.cache_data
# def load_data():
#     try:
#         df_loaded = pd.read_csv("https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/peptides/peptide_groups.csv")
#         if 'group' not in df_loaded.columns or 'peptide' not in df_loaded.columns:
#             st.error("Error: DataFrame must contain 'group' and 'peptide' columns.")
#             return pd.DataFrame() # Return empty DataFrame on error
#         # Clean up peptides: strip whitespace and drop any empty ones
#         df_loaded['peptide'] = df_loaded['peptide'].astype(str).str.strip()
#         df_loaded = df_loaded[df_loaded['peptide'] != ''].copy()
        
#         # Ensure 'group' column is also clean and not empty after dropping NaNs
#         df_loaded['group'] = df_loaded['group'].astype(str).str.strip()
#         df_loaded = df_loaded[df_loaded['group'] != ''].copy()

#         return df_loaded
#     except Exception as e:
#         st.error(f"Error loading DataFrame from URL: {e}")
#         # Fallback to dummy data for demonstration if URL fails
#         data = {
#             'group': ['group_A']*50 + ['group_B']*50 + ['group_C']*50,
#             'peptide': [
#                 'AYLDVFNKF', 'VYIDVFNWD', 'AALDVFNKF', 'PYLDVFNKW', 'AYLDFNKKF',
#             ]*10 + [
#                 'VYIDVFNWD', 'GYIDVFNWD', 'TYLDVFNWS', 'VYIDVENWA', 'AYLDVFNKE',
#             ]*10 + [
#                 'PPLDVFNKA', 'QYLDVFNKC', 'RYIDVFNKD', 'AYLDVFNKG', 'AYLDVFNKH'
#             ]*10
#         }
#         import random
#         import string
#         def generate_random_peptide(length=9):
#             return ''.join(random.choice(string.ascii_uppercase) for _ in range(length))
#         for i in range(len(data['peptide'])):
#             if random.random() < 0.2: # Introduce some randomness
#                 data['peptide'][i] = generate_random_peptide()
#         st.warning("Using a dummy DataFrame for demonstration (real CSV failed to load).")
#         return pd.DataFrame(data)

# # --- Streamlit App Layout ---
# st.set_page_config(layout="wide")
# st.title("Peptide Group Similarity Dashboard")
# st.write("This dashboard allows you to explore peptide group similarities based on sequence matching.")

# # --- Data Source Selection ---
# data_source_option = st.radio(
#     "Choose data source:",
#     ("Use existing 'Peptide Groups' data", "Upload a custom CSV file"),
#     key="data_source_radio"
# )

# df = None
# if data_source_option == "Use existing 'Peptide Groups' data":
#     df = load_data()
#     if not df.empty:
#         st.success("Default 'Peptide Groups' data loaded successfully!")
#     else:
#         st.error("Could not load default peptide data. Please try uploading a custom CSV.")

# elif data_source_option == "Upload a custom CSV file":
#     uploaded_file = st.file_uploader("Upload your own CSV file", type="csv", key="file_uploader")
#     if uploaded_file is not None:
#         try:
#             df = pd.read_csv(uploaded_file)
#             # Basic validation for peptide data if custom CSV is uploaded
#             if 'group' not in df.columns or 'peptide' not in df.columns:
#                 st.warning("Warning: Uploaded CSV does not contain 'group' and 'peptide' columns. Some features may not work as expected.")
            
#             # Clean up peptides for uploaded data too
#             df['peptide'] = df['peptide'].astype(str).str.strip()
#             df = df[df['peptide'] != ''].copy()
#             # Clean up groups for uploaded data too
#             df['group'] = df['group'].astype(str).str.strip()
#             df = df[df['group'] != ''].copy()

#             st.success("Custom CSV file uploaded and processed!")

#         except Exception as e:
#             st.error(f"Error reading uploaded file: {e}. Please ensure it's a valid CSV.")
#             df = None
#     else:
#         st.info("Please upload a CSV file to proceed.")


# # --- Main Dashboard Logic (Conditional on Data Loading) ---
# if df is not None and not df.empty:
#     st.sidebar.header("Data Overview")
#     st.sidebar.write(f"**Total Peptides:** {len(df)}")
#     st.sidebar.write(f"**Number of Groups:** {df['group'].nunique()}")
    
#     st.subheader("Data Preview")
#     st.dataframe(df.head())

#     # --- Section 1: Analyze Reference Peptide Similarity ---
#     st.header("1. Reference Peptide Similarity Analysis")
#     st.write("Determine how similar a chosen peptide is to all peptides within each group.")

#     st.sidebar.subheader("Reference Peptide Configuration")

#     # --- NEW: Choose source for reference peptide ---
#     ref_peptide_source_option = st.sidebar.radio(
#         "How to specify Reference Peptide?",
#         ("Select from existing peptides", "Enter custom peptide"),
#         key="ref_peptide_source_radio"
#     )

#     selected_reference_peptide = None # Initialize to None

#     # Get the unique characters (alphabet) from all peptides in the loaded data
#     all_peptides_str = "".join(df['peptide'].dropna().tolist())
#     peptide_alphabet = sorted(list(set(char.upper() for char in all_peptides_str if char.isalpha())))
    
#     # Display allowed characters in sidebar
#     if peptide_alphabet:
#         st.sidebar.info(f"**Allowed peptide characters:** {' '.join(peptide_alphabet)}")
#     else:
#         st.sidebar.warning("Could not determine a peptide alphabet from loaded data.")


#     if ref_peptide_source_option == "Select from existing peptides":
#         example_peptides = df['peptide'].dropna().unique().tolist()
#         if not example_peptides:
#             st.sidebar.warning("No peptides found in data to select from.")
#             selected_reference_peptide = None
#         else:
#             # Try to set a common default, otherwise pick the first
#             default_ref_peptide_val = 'AYLDVFNKF'
#             default_index = example_peptides.index(default_ref_peptide_val) if default_ref_peptide_val in example_peptides else 0

#             selected_reference_peptide = st.sidebar.selectbox(
#                 "Choose from existing peptides:",
#                 options=example_peptides,
#                 index=default_index,
#                 key="select_existing_ref_peptide"
#             )

#     elif ref_peptide_source_option == "Enter custom peptide":
#         # Determine a reasonable max length for the input field from existing peptides
#         max_peptide_len = df['peptide'].str.len().max() if not df.empty else 20 # Default if df is empty

#         custom_peptide_input = st.sidebar.text_input(
#             "Enter your custom peptide (e.g., 'AYLDVFNKF'):",
#             value='AYLDVFNKF', # Default value for custom input
#             max_chars=max_peptide_len,
#             key="custom_ref_peptide_input"
#         )
#         custom_peptide_input_clean = custom_peptide_input.strip().upper() # Clean and uppercase for validation

#         # --- Validation Logic for Custom Peptide ---
#         if custom_peptide_input_clean:
#             # 1. Check if it's purely alphabetical
#             if not custom_peptide_input_clean.isalpha():
#                 st.sidebar.error("Custom peptide must contain only alphabetical characters (A-Z, a-z).")
#                 selected_reference_peptide = None
#             else:
#                 # 2. Check if all characters are in the derived alphabet
#                 is_valid_chars = True
#                 invalid_chars = []
#                 for char in custom_peptide_input_clean:
#                     if char not in peptide_alphabet:
#                         is_valid_chars = False
#                         invalid_chars.append(char)
                
#                 if not is_valid_chars:
#                     st.sidebar.error(f"Invalid characters detected: {' '.join(set(invalid_chars))}. Only these are allowed: {' '.join(peptide_alphabet)}")
#                     selected_reference_peptide = None # Invalidate peptide if not valid
#                 else:
#                     selected_reference_peptide = custom_peptide_input_clean # Valid peptide, use uppercase for consistency
#         else:
#             st.sidebar.info("Please enter a peptide sequence.")
#             selected_reference_peptide = None # No peptide entered

#     # Filter groups for analysis
#     all_groups = df['group'].dropna().unique().tolist()
#     if not all_groups: # Handle case of no groups
#         st.warning("No peptide groups found in the loaded data to analyze.")
#     elif selected_reference_peptide: # Only proceed if a valid reference peptide is chosen
#         selected_groups_for_ref = st.multiselect(
#             "Select Groups to Analyze (for reference peptide):",
#             options=all_groups,
#             default=all_groups, # Default to all groups
#             key="ref_peptide_group_multiselect"
#         )

#         if selected_groups_for_ref:
#             # Filter df for selected groups before calculating scores
#             df_filtered_for_ref = df[df['group'].isin(selected_groups_for_ref)]
            
#             # Calculate scores
#             df_with_scores = calculate_similarity_scores(selected_reference_peptide, df_filtered_for_ref)
#             valid_scores_df = df_with_scores.dropna(subset=['similarity_score'])

#             if not valid_scores_df.empty:
#                 # Box Plot
#                 st.subheader(f"Similarity Score Distribution to '{selected_reference_peptide}'")
#                 fig1, ax1 = plt.subplots(figsize=(12, 7))
#                 sns.boxplot(x='group', y='similarity_score', data=valid_scores_df, ax=ax1, palette='viridis')
#                 ax1.set_title(f'Similarity Score Distribution to "{selected_reference_peptide}" Across Selected Groups')
#                 ax1.set_xlabel('Peptide Group')
#                 ax1.set_ylabel(f'Similarity Score to "{selected_reference_peptide}"')
#                 ax1.set_ylim(0, 1) # Scores are 0 to 1
#                 ax1.grid(axis='y', linestyle='--', alpha=0.7)
#                 st.pyplot(fig1)

#                 # Quantitative Spread Stats
#                 st.subheader("Quantitative Spread within Each Group:")
#                 group_spread_stats = valid_scores_df.groupby('group')['similarity_score'].agg(
#                     mean_score='mean',
#                     median_score='median',
#                     std_dev_score='std',
#                     iqr_score=lambda x: x.quantile(0.75) - x.quantile(0.25),
#                     count='count'
#                 ).dropna()
#                 st.dataframe(group_spread_stats.style.format("{:.4f}"))
#             else:
#                 st.info(f"No valid similarity scores found for '{selected_reference_peptide}' within the selected groups. This might be due to peptide length mismatch or group content.")
#         else:
#             st.info("Please select at least one group to analyze against the reference peptide.")
#     else:
#         st.info("Please select or enter a valid reference peptide to proceed with analysis.")


#     # --- Section 2: Explore Inter-Group Similarity ---
#     st.header("2. Inter-Group Similarity Analysis")
#     st.write("Identify and visualize how similar different peptide groups are to each other.")

#     # Calculate all group-to-group similarities
#     group_similarity_results = determine_most_similar_groups(df.copy())

#     if group_similarity_results['group_pair'] and group_similarity_results['all_pair_scores']:
#         # Sort for clear presentation
#         sorted_pair_scores_list = sorted(group_similarity_results['all_pair_scores'].items(), key=lambda item: item[1], reverse=True)
#         pair_scores_df = pd.DataFrame(
#             [(f"{pair[0]} vs {pair[1]}", score) for pair, score in sorted_pair_scores_list],
#             columns=['Group Pair', 'Average Similarity Score']
#         )

#         st.subheader("Average Similarity Scores Between All Group Pairs:")
#         num_top_pairs = st.slider("Show Top N Most Similar Pairs:", 1, len(pair_scores_df), min(5, len(pair_scores_df)), key="num_top_pairs_slider")
        
#         # Display top N pairs in a bar chart
#         fig2, ax2 = plt.subplots(figsize=(12, max(6, num_top_pairs * 0.8))) # Adjust height dynamically
#         sns.barplot(x='Average Similarity Score', y='Group Pair', data=pair_scores_df.head(num_top_pairs), palette='coolwarm', ax=ax2)
#         ax2.set_title(f'Top {num_top_pairs} Most Similar Peptide Group Pairs')
#         ax2.set_xlabel('Average Similarity Score (0-1)')
#         ax2.set_ylabel('Peptide Group Pair')
#         ax2.set_xlim(0, 1)
#         ax2.grid(axis='x', linestyle='--', alpha=0.7)
#         st.pyplot(fig2)

#         st.success(f"**The single most similar pair is: {group_similarity_results['group_pair'][0]} and {group_similarity_results['group_pair'][1]}** (Average Score: {group_similarity_results['average_similarity_score']:.4f})")

#         # Detailed Heatmap for Selected Pair
#         st.subheader("Detailed Pairwise Similarity Heatmap")
        
#         # Create a list of available group pairs for the selectbox
#         available_pairs_for_heatmap = [f"{pair[0]} vs {pair[1]}" for pair, score in sorted_pair_scores_list]
#         selected_heatmap_pair_str = st.selectbox(
#             "Select a Group Pair for Detailed Heatmap:",
#             options=available_pairs_for_heatmap,
#             index=0, # Default to the most similar pair
#             key="heatmap_pair_selectbox"
#         )

#         # Parse the selected string back into group names
#         selected_group1_name, selected_group2_name = selected_heatmap_pair_str.split(' vs ')

#         peptides_g1 = df[df['group'] == selected_group1_name]['peptide'].dropna().tolist()
#         peptides_g2 = df[df['group'] == selected_group2_name]['peptide'].dropna().tolist()

#         # Slider to limit the number of peptides shown on the heatmap
#         max_heatmap_peptides = st.slider(
#             "Limit # of peptides for heatmap (adjust for clarity):",
#             min_value=5, max_value=min(50, max(len(peptides_g1), len(peptides_g2), 10)), # Max up to 50 or max group size
#             value=min(20, max(len(peptides_g1), len(peptides_g2), 10)), # Default to 20 or group size
#             key="heatmap_peptide_limit"
#         )
#         peptides_g1_sampled = peptides_g1[:max_heatmap_peptides]
#         peptides_g2_sampled = peptides_g2[:max_heatmap_peptides]

#         similarity_matrix_inter_group = pd.DataFrame(
#             np.nan,
#             index=peptides_g1_sampled,
#             columns=peptides_g2_sampled,
#             dtype=float
#         )

#         # Populate heatmap matrix
#         for p1 in peptides_g1_sampled:
#             for p2 in peptides_g2_sampled:
#                 if not isinstance(p1, str) or not p1.strip() or not isinstance(p2, str) or not p2.strip():
#                     similarity_matrix_inter_group.loc[p1, p2] = np.nan
#                     continue
#                 try:
#                     score = calculate_peptide_similarity(p1, p2)
#                     similarity_matrix_inter_group.loc[p1, p2] = score
#                 except (ValueError, TypeError):
#                     similarity_matrix_inter_group.loc[p1, p2] = np.nan
#                     continue

#         # Check if the matrix is completely empty (all NaNs) before plotting
#         if not similarity_matrix_inter_group.empty and not similarity_matrix_inter_group.isna().all().all():
#             fig4, ax4 = plt.subplots(figsize=(max(10, len(peptides_g2_sampled)*0.8), max(8, len(peptides_g1_sampled)*0.8)))
#             sns.heatmap(
#                 similarity_matrix_inter_group,
#                 cmap='viridis',
#                 annot=True if (len(peptides_g1_sampled) < 15 and len(peptides_g2_sampled) < 15) else False,
#                 fmt=".2f",
#                 cbar_kws={'label': 'Similarity Score (0-1)'},
#                 ax=ax4
#             )
#             ax4.set_title(f'Pairwise Similarity Between {selected_group1_name} and {selected_group2_name}')
#             ax4.set_xlabel(f'Peptides in {selected_group2_name}')
#             ax4.set_ylabel(f'Peptides in {selected_group1_name}')
#             plt.tight_layout() # Adjust layout to prevent labels from overlapping
#             st.pyplot(fig4)
#         else:
#             st.info("No valid data or peptides for this heatmap. Adjust selections or check data content (e.g., peptide strings might be empty).")
#     else:
#         st.warning("Not enough distinct peptide groups to perform inter-group similarity analysis (need at least 2), or no valid peptide data loaded.")

# else:
#     st.info("Please load peptide data using the options above to enable analysis.")










# import streamlit as st
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import itertools
# import re # Import regex for parsing mutant string

# # --- Peptide Similarity Calculation Function (re-use from previous code) ---
# def calculate_peptide_similarity(peptide_1: str, peptide_2: str) -> float:
#     # Ensure peptides are processed as strings and stripped
#     peptide_1 = str(peptide_1).strip()
#     peptide_2 = str(peptide_2).strip()

#     if not peptide_1 or not peptide_2: # Handle empty strings
#         return 0.0

#     # It's a design choice to return 0.0 for different lengths
#     # or raise an error if strict equality is needed.
#     # For a dashboard, 0.0 is often more user-friendly.
#     if len(peptide_1) != len(peptide_2):
#         return 0.0
            
#     num_matched = 0
#     for i in range(len(peptide_1)):
#         if peptide_1[i] == peptide_2[i]:
#             num_matched += 1
            
#     return num_matched / len(peptide_1)

# # --- Function to calculate similarity scores for a reference peptide (re-use) ---
# @st.cache_data # Cache this function to speed up reruns
# def calculate_similarity_scores(reference_peptide: str, df: pd.DataFrame) -> pd.DataFrame:
#     df_copy = df.copy()
#     df_copy['similarity_score'] = float('nan')
    
#     # Ensure reference_peptide is clean
#     reference_peptide = str(reference_peptide).strip()
#     if not reference_peptide: # If reference peptide is empty, no scores can be calculated
#         return df_copy

#     for index, row in df_copy.iterrows():
#         peptide_x = str(row['peptide']).strip() # Ensure it's a string and strip whitespace
#         if not peptide_x: # Skip empty peptide strings in the dataframe
#             continue
#         try:
#             score = calculate_peptide_similarity(reference_peptide, peptide_x)
#             df_copy.at[index, 'similarity_score'] = score
#         except (ValueError, TypeError): # Catch errors from calculate_peptide_similarity
#             continue
#     return df_copy

# # --- Function to determine most similar groups (re-use) ---
# @st.cache_data # Cache this function as it can be computationally intensive
# def determine_most_similar_groups(df: pd.DataFrame) -> dict:
#     unique_groups = df['group'].dropna().unique().tolist()
#     if len(unique_groups) < 2: # Need at least two groups to find a pair
#         return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}

#     group_pairs = list(itertools.combinations(unique_groups, 2))
#     all_group_pair_scores = {}

#     for group1_name, group2_name in group_pairs:
#         peptides_group1 = df[df['group'] == group1_name]['peptide'].dropna().tolist()
#         peptides_group2 = df[df['group'] == group2_name]['peptide'].dropna().tolist()
#         pairwise_scores_for_current_pair = []

#         if not peptides_group1 or not peptides_group2: # Skip if a group is empty
#             all_group_pair_scores[(group1_name, group2_name)] = 0.0
#             continue

#         for p1 in peptides_group1:
#             if not isinstance(p1, str) or not p1.strip(): continue # Skip non-strings or empty
#             for p2 in peptides_group2:
#                 if not isinstance(p2, str) or not p2.strip(): continue # Skip non-strings or empty
#                 try:
#                     score = calculate_peptide_similarity(p1.strip(), p2.strip())
#                     pairwise_scores_for_current_pair.append(score)
#                 except (ValueError, TypeError): # Catch errors from calculate_peptide_similarity
#                     continue
        
#         if pairwise_scores_for_current_pair:
#             avg_score = sum(pairwise_scores_for_current_pair) / len(pairwise_scores_for_current_pair)
#             all_group_pair_scores[(group1_name, group2_name)] = avg_score
#         else:
#             all_group_pair_scores[(group1_name, group2_name)] = 0.0

#     if not all_group_pair_scores: # If no scores were calculated at all
#         return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}
        
#     most_similar_pair = max(all_group_pair_scores.items(), key=lambda item: item[1])

#     return {
#         'group_pair': most_similar_pair[0],
#         'average_similarity_score': most_similar_pair[1],
#         'all_pair_scores': all_group_pair_scores
#     }

# # --- Data Loading Functions ---
# @st.cache_data
# def load_peptide_groups_data():
#     try:
#         df_loaded = pd.read_csv("https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/peptides/peptide_groups.csv")
#         if 'group' not in df_loaded.columns or 'peptide' not in df_loaded.columns:
#             st.error("Error: 'peptide_groups.csv' must contain 'group' and 'peptide' columns.")
#             return pd.DataFrame() # Return empty DataFrame on error
#         # Clean up peptides: strip whitespace and drop any empty ones
#         df_loaded['peptide'] = df_loaded['peptide'].astype(str).str.strip()
#         df_loaded = df_loaded[df_loaded['peptide'] != ''].copy()
        
#         # Ensure 'group' column is also clean and not empty after dropping NaNs
#         df_loaded['group'] = df_loaded['group'].astype(str).str.strip()
#         df_loaded = df_loaded[df_loaded['group'] != ''].copy()

#         return df_loaded
#     except Exception as e:
#         st.error(f"Error loading default 'peptide_groups.csv' from URL: {e}")
#         st.warning("Using a dummy 'Peptide Groups' DataFrame for demonstration due to loading failure.")
#         # Fallback to dummy data for demonstration if URL fails
#         data = {
#             'group': ['group_A']*50 + ['group_B']*50 + ['group_C']*50,
#             'peptide': [
#                 'AYLDVFNKF', 'VYIDVFNWD', 'AALDVFNKF', 'PYLDVFNKW', 'AYLDFNKKF',
#             ]*10 + [
#                 'VYIDVFNWD', 'GYIDVFNWD', 'TYLDVFNWS', 'VYIDVENWA', 'AYLDVFNKE',
#             ]*10 + [
#                 'PPLDVFNKA', 'QYLDVFNKC', 'RYIDVFNKD', 'AYLDVFNKG', 'AYLDVFNKH'
#             ]*10
#         }
#         import random
#         import string
#         def generate_random_peptide(length=9):
#             return ''.join(random.choice(string.ascii_uppercase) for _ in range(length))
#         for i in range(len(data['peptide'])):
#             if random.random() < 0.2: # Introduce some randomness
#                 data['peptide'][i] = generate_random_peptide()
#         return pd.DataFrame(data)

# @st.cache_data
# def load_dms_data():
#     try:
#         dms_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/dms.csv')
#         if 'mutant' not in dms_df.columns or 'DMS_score' not in dms_df.columns:
#             st.error("Error: 'dms.csv' must contain 'mutant' and 'DMS_score' columns.")
#             return pd.DataFrame()
        
#         # Parse the 'mutant' column: e.g., A123G -> Original: A, Position: 123, New: G
#         mutant_pattern = re.compile(r'([A-Z])(\d+)([A-Z])')
        
#         # Apply regex and handle potential failures gracefully
#         parsed_results = dms_df['mutant'].apply(lambda x: mutant_pattern.match(str(x)))
        
#         # Extract components; use .get(idx, None) for safety if match fails
#         dms_df['original_aa'] = [m.group(1) if m else None for m in parsed_results]
#         dms_df['position'] = [int(m.group(2)) if m and m.group(2).isdigit() else None for m in parsed_results]
#         dms_df['new_aa'] = [m.group(3) if m else None for m in parsed_results]
        
#         # Drop rows where parsing failed or essential columns are missing
#         dms_df.dropna(subset=['original_aa', 'position', 'new_aa', 'DMS_score'], inplace=True)
#         dms_df['position'] = dms_df['position'].astype(int) # Ensure position is integer
        
#         return dms_df
#     except Exception as e:
#         st.error(f"Error loading 'dms.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
#         # Do NOT return an empty DataFrame immediately here,
#         # let the tab_dms content handle the empty DataFrame check.
#         # This allows the rest of the tab to render its messages.
#         return pd.DataFrame() # Still return empty DataFrame to signify no data loaded


# # --- Streamlit App Layout ---
# st.set_page_config(layout="wide")
# st.title("Protein Analysis Dashboard")
# st.write("Explore peptide group similarities and deep mutational scanning (DMS) data.")

# # --- Tabbed Interface for different analyses ---
# # Explicitly set a default tab if you want one to be active on load
# tab_peptide, tab_dms = st.tabs(["Peptide Group Analysis", "Deep Mutational Scanning (DMS) Analysis"])

# # --- Peptide Group Analysis Tab ---
# with tab_peptide:
#     st.header("Peptide Group Similarity Dashboard")
#     st.write("This dashboard allows you to explore peptide group similarities based on sequence matching.")

#     # --- Data Source Selection for Peptide Analysis ---
#     data_source_option = st.radio(
#         "Choose data source for Peptide Analysis:",
#         ("Use existing 'Peptide Groups' data", "Upload a custom CSV file"),
#         key="peptide_data_source_radio"
#     )

#     df_peptide = None # Use a distinct name for peptide DataFrame
#     if data_source_option == "Use existing 'Peptide Groups' data":
#         df_peptide = load_peptide_groups_data()
#         if not df_peptide.empty:
#             st.success("Default 'Peptide Groups' data loaded successfully!")
#         # else: st.error message is handled by load_peptide_groups_data itself
            
#     elif data_source_option == "Upload a custom CSV file":
#         uploaded_file = st.file_uploader("Upload your own CSV file for Peptide Analysis", type="csv", key="peptide_file_uploader")
#         if uploaded_file is not None:
#             try:
#                 df_peptide = pd.read_csv(uploaded_file)
#                 if 'group' not in df_peptide.columns or 'peptide' not in df_peptide.columns:
#                     st.warning("Warning: Uploaded CSV does not contain 'group' and 'peptide' columns. Some features may not work.")
                
#                 df_peptide['peptide'] = df_peptide['peptide'].astype(str).str.strip()
#                 df_peptide = df_peptide[df_peptide['peptide'] != ''].copy()
#                 df_peptide['group'] = df_peptide['group'].astype(str).str.strip()
#                 df_peptide = df_peptide[df_peptide['group'] != ''].copy()

#                 st.success("Custom CSV file uploaded and processed for Peptide Analysis!")

#             except Exception as e:
#                 st.error(f"Error reading uploaded file: {e}. Please ensure it's a valid CSV.")
#                 df_peptide = None
#         else:
#             st.info("Please upload a CSV file to proceed with Peptide Analysis.")


#     # --- Main Peptide Dashboard Logic (Conditional on Data Loading) ---
#     if df_peptide is not None and not df_peptide.empty:
#         st.sidebar.header("Peptide Data Overview")
#         st.sidebar.write(f"**Total Peptides:** {len(df_peptide)}")
#         st.sidebar.write(f"**Number of Groups:** {df_peptide['group'].nunique()}")
        
#         st.subheader("Peptide Data Preview")
#         st.dataframe(df_peptide.head())

#         # --- Section 1: Analyze Reference Peptide Similarity ---
#         st.header("1. Reference Peptide Similarity Analysis")
#         st.write("Determine how similar a chosen peptide is to all peptides within each group.")

#         st.sidebar.subheader("Reference Peptide Configuration")

#         ref_peptide_source_option = st.sidebar.radio(
#             "How to specify Reference Peptide?",
#             ("Select from existing peptides", "Enter custom peptide"),
#             key="ref_peptide_source_radio"
#         )

#         selected_reference_peptide = None 

#         all_peptides_str = "".join(df_peptide['peptide'].dropna().tolist())
#         peptide_alphabet = sorted(list(set(char.upper() for char in all_peptides_str if char.isalpha())))
        
#         if peptide_alphabet:
#             st.sidebar.info(f"**Allowed peptide characters:** {' '.join(peptide_alphabet)}")
#         else:
#             st.sidebar.warning("Could not determine a peptide alphabet from loaded peptide data. Please ensure 'peptide' column contains valid sequences.")


#         if ref_peptide_source_option == "Select from existing peptides":
#             example_peptides = df_peptide['peptide'].dropna().unique().tolist()
#             if not example_peptides:
#                 st.sidebar.warning("No peptides found in data to select from.")
#                 selected_reference_peptide = None
#             else:
#                 default_ref_peptide_val = 'AYLDVFNKF'
#                 # Ensure the default_ref_peptide_val is in example_peptides before finding index
#                 default_index = example_peptides.index(default_ref_peptide_val) if default_ref_peptide_val in example_peptides else 0
                
#                 selected_reference_peptide = st.sidebar.selectbox(
#                     "Choose from existing peptides:",
#                     options=example_peptides,
#                     index=default_index, # Fixed to use the correct default_index
#                     key="select_existing_ref_peptide"
#                 )

#         elif ref_peptide_source_option == "Enter custom peptide":
#             max_peptide_len = df_peptide['peptide'].str.len().max() if not df_peptide.empty else 20

#             custom_peptide_input = st.sidebar.text_input(
#                 "Enter your custom peptide (e.g., 'AYLDVFNKF'):",
#                 value='AYLDVFNKF', # Default value for custom input
#                 max_chars=max_peptide_len,
#                 key="custom_ref_peptide_input"
#             )
#             custom_peptide_input_clean = custom_peptide_input.strip().upper()

#             if custom_peptide_input_clean:
#                 if not custom_peptide_input_clean.isalpha():
#                     st.sidebar.error("Custom peptide must contain only alphabetical characters (A-Z, a-z).")
#                     selected_reference_peptide = None
#                 else:
#                     is_valid_chars = True
#                     invalid_chars = []
#                     for char in custom_peptide_input_clean:
#                         if char not in peptide_alphabet:
#                             is_valid_chars = False
#                             invalid_chars.append(char)
                    
#                     if not is_valid_chars:
#                         st.sidebar.error(f"Invalid characters detected: {' '.join(set(invalid_chars))}. Only these are allowed: {' '.join(peptide_alphabet)}")
#                         selected_reference_peptide = None
#                     else:
#                         selected_reference_peptide = custom_peptide_input_clean
#             else:
#                 st.sidebar.info("Please enter a peptide sequence.")
#                 selected_reference_peptide = None

#         all_groups = df_peptide['group'].dropna().unique().tolist()
#         if not all_groups:
#             st.warning("No peptide groups found in the loaded data to analyze.")
#         elif selected_reference_peptide: # Only proceed if a valid reference peptide is chosen
#             selected_groups_for_ref = st.multiselect(
#                 "Select Groups to Analyze (for reference peptide):",
#                 options=all_groups,
#                 default=all_groups,
#                 key="ref_peptide_group_multiselect"
#             )

#             if selected_groups_for_ref:
#                 df_filtered_for_ref = df_peptide[df_peptide['group'].isin(selected_groups_for_ref)]
                
#                 df_with_scores = calculate_similarity_scores(selected_reference_peptide, df_filtered_for_ref)
#                 valid_scores_df = df_with_scores.dropna(subset=['similarity_score'])

#                 if not valid_scores_df.empty:
#                     st.subheader(f"Similarity Score Distribution to '{selected_reference_peptide}'")
#                     fig1, ax1 = plt.subplots(figsize=(12, 7))
#                     sns.boxplot(x='group', y='similarity_score', data=valid_scores_df, ax=ax1, palette='viridis')
#                     ax1.set_title(f'Similarity Score Distribution to "{selected_reference_peptide}" Across Selected Groups')
#                     ax1.set_xlabel('Peptide Group')
#                     ax1.set_ylabel(f'Similarity Score to "{selected_reference_peptide}"')
#                     ax1.set_ylim(0, 1)
#                     ax1.grid(axis='y', linestyle='--', alpha=0.7)
#                     st.pyplot(fig1)

#                     st.subheader("Quantitative Spread within Each Group:")
#                     group_spread_stats = valid_scores_df.groupby('group')['similarity_score'].agg(
#                         mean_score='mean',
#                         median_score='median',
#                         std_dev_score='std',
#                         iqr_score=lambda x: x.quantile(0.75) - x.quantile(0.25),
#                         count='count'
#                     ).dropna()
#                     st.dataframe(group_spread_stats.style.format("{:.4f}"))
#                 else:
#                     st.info(f"No valid similarity scores found for '{selected_reference_peptide}' within the selected groups. This might be due to peptide length mismatch or group content.")
#             else:
#                 st.info("Please select at least one group to analyze against the reference peptide.")
#         else:
#             st.info("Please select or enter a valid reference peptide to proceed with analysis.")

#         # --- Section 2: Explore Inter-Group Similarity ---
#         st.header("2. Inter-Group Similarity Analysis")
#         st.write("Identify and visualize how similar different peptide groups are to each other.")

#         group_similarity_results = determine_most_similar_groups(df_peptide.copy())

#         if group_similarity_results['group_pair'] and group_similarity_results['all_pair_scores']:
#             sorted_pair_scores_list = sorted(group_similarity_results['all_pair_scores'].items(), key=lambda item: item[1], reverse=True)
#             pair_scores_df = pd.DataFrame(
#                 [(f"{pair[0]} vs {pair[1]}", score) for pair, score in sorted_pair_scores_list],
#                 columns=['Group Pair', 'Average Similarity Score']
#             )

#             st.subheader("Average Similarity Scores Between All Group Pairs:")
#             num_top_pairs = st.slider("Show Top N Most Similar Pairs:", 1, len(pair_scores_df), min(5, len(pair_scores_df)), key="num_top_pairs_slider")
            
#             fig2, ax2 = plt.subplots(figsize=(12, max(6, num_top_pairs * 0.8)))
#             sns.barplot(x='Average Similarity Score', y='Group Pair', data=pair_scores_df.head(num_top_pairs), palette='coolwarm', ax=ax2)
#             ax2.set_title(f'Top {num_top_pairs} Most Similar Peptide Group Pairs')
#             ax2.set_xlabel('Average Similarity Score (0-1)')
#             ax2.set_ylabel('Peptide Group Pair')
#             ax2.set_xlim(0, 1)
#             ax2.grid(axis='x', linestyle='--', alpha=0.7)
#             st.pyplot(fig2)

#             st.success(f"**The single most similar pair is: {group_similarity_results['group_pair'][0]} and {group_similarity_results['group_pair'][1]}** (Average Score: {group_similarity_results['average_similarity_score']:.4f})")

#             st.subheader("Detailed Pairwise Similarity Heatmap")
            
#             available_pairs_for_heatmap = [f"{pair[0]} vs {pair[1]}" for pair, score in sorted_pair_scores_list]
#             selected_heatmap_pair_str = st.selectbox(
#                 "Select a Group Pair for Detailed Heatmap:",
#                 options=available_pairs_for_heatmap,
#                 index=0,
#                 key="heatmap_pair_selectbox"
#             )

#             selected_group1_name, selected_group2_name = selected_heatmap_pair_str.split(' vs ')

#             peptides_g1 = df_peptide[df_peptide['group'] == selected_group1_name]['peptide'].dropna().tolist()
#             peptides_g2 = df_peptide[df_peptide['group'] == selected_group2_name]['peptide'].dropna().tolist()

#             max_heatmap_peptides = st.slider(
#                 "Limit # of peptides for heatmap (adjust for clarity):",
#                 min_value=5, max_value=min(50, max(len(peptides_g1), len(peptides_g2), 10)),
#                 value=min(20, max(len(peptides_g1), len(peptides_g2), 10)),
#                 key="heatmap_peptide_limit"
#             )
#             peptides_g1_sampled = peptides_g1[:max_heatmap_peptides]
#             peptides_g2_sampled = peptides_g2[:max_heatmap_peptides]

#             similarity_matrix_inter_group = pd.DataFrame(
#                 np.nan,
#                 index=peptides_g1_sampled,
#                 columns=peptides_g2_sampled,
#                 dtype=float
#             )

#             for p1 in peptides_g1_sampled:
#                 for p2 in peptides_g2_sampled:
#                     if not isinstance(p1, str) or not p1.strip() or not isinstance(p2, str) or not p2.strip():
#                         similarity_matrix_inter_group.loc[p1, p2] = np.nan
#                         continue
#                     try:
#                         score = calculate_peptide_similarity(p1, p2)
#                         similarity_matrix_inter_group.loc[p1, p2] = score
#                     except (ValueError, TypeError):
#                         similarity_matrix_inter_group.loc[p1, p2] = np.nan
#                         continue

#             if not similarity_matrix_inter_group.empty and not similarity_matrix_inter_group.isna().all().all():
#                 fig4, ax4 = plt.subplots(figsize=(max(10, len(peptides_g2_sampled)*0.8), max(8, len(peptides_g1_sampled)*0.8)))
#                 sns.heatmap(
#                     similarity_matrix_inter_group,
#                     cmap='viridis',
#                     annot=True if (len(peptides_g1_sampled) < 15 and len(peptides_g2_sampled) < 15) else False,
#                     fmt=".2f",
#                     cbar_kws={'label': 'Similarity Score (0-1)'},
#                     ax=ax4
#                 )
#                 ax4.set_title(f'Pairwise Similarity Between {selected_group1_name} and {selected_group2_name}')
#                 ax4.set_xlabel(f'Peptides in {selected_group2_name}')
#                 ax4.set_ylabel(f'Peptides in {selected_group1_name}')
#                 plt.tight_layout()
#                 st.pyplot(fig4)
#             else:
#                 st.info("No valid data or peptides for this heatmap. Adjust selections or check data content (e.g., peptide strings might be empty).")
#         else:
#             st.warning("Not enough distinct peptide groups to perform inter-group similarity analysis (need at least 2), or no valid peptide data loaded.")

#     else:
#         st.info("Please load peptide data using the options above to enable analysis.")


# # --- Deep Mutational Scanning (DMS) Analysis Tab ---
# with tab_dms:
#     st.header("Deep Mutational Scanning (DMS) Analysis")
#     st.write("Visualize DMS scores as a heatmap to identify positions and mutations affecting protein fitness/function.")

#     dms_df = load_dms_data()

#     if not dms_df.empty:
#         st.success("DMS data loaded successfully!")
#         st.subheader("DMS Data Preview")
#         st.dataframe(dms_df.head())

#         st.subheader("DMS Heatmap Configuration")

#         # Get unique positions and new amino acids
#         all_positions = sorted(dms_df['position'].unique())
#         all_new_aas = sorted(dms_df['new_aa'].unique())

#         if not all_positions or not all_new_aas:
#             st.warning("Not enough data to generate heatmap (missing positions or new amino acids after parsing).")
#         else:
#             # Position Range Slider
#             min_pos, max_pos = min(all_positions), max(all_positions)
#             selected_pos_range = st.slider(
#                 "Select Position Range:",
#                 min_value=min_pos,
#                 max_value=max_pos,
#                 value=(min_pos, max_pos),
#                 key="dms_pos_range_slider"
#             )

#             # Filter data based on selected position range
#             filtered_dms_df = dms_df[(dms_df['position'] >= selected_pos_range[0]) & 
#                                      (dms_df['position'] <= selected_pos_range[1])]

#             # Aggregate scores for the heatmap
#             heatmap_data = filtered_dms_df.pivot_table(
#                 values='DMS_score',
#                 index='new_aa',
#                 columns='position',
#                 aggfunc='mean'
#             ).fillna(np.nan) # Use NaN for missing values in heatmap to distinguish from 0 scores

#             # Reindex to ensure all unique new_aas and all positions in range are included
#             # This ensures consistent heatmap structure even if no data for certain positions/aas
#             heatmap_data = heatmap_data.reindex(index=all_new_aas, columns=range(min_pos, max_pos + 1), fill_value=np.nan)
            
#             # --- Heatmap Display ---
#             st.subheader("DMS Score Heatmap: Average Score per Mutation")
#             st.write("Colors represent the average DMS score for mutations. Brighter (yellow) indicates higher scores, suggesting mutations are more tolerated at that position to that new amino acid.")

#             # Color scale control (optional, but good for interpretation)
#             # Calculate initial vmin/vmax from *actual* data, not NaNs
#             valid_scores_flat = heatmap_data.stack().dropna() # Get all valid scores
#             if not valid_scores_flat.empty:
#                 default_vmin = valid_scores_flat.min()
#                 default_vmax = valid_scores_flat.max()
#             else:
#                 default_vmin = 0.0
#                 default_vmax = 1.0 # Fallback default if no valid scores

#             col1, col2 = st.columns(2)
#             with col1:
#                 heatmap_vmin = st.number_input("Heatmap Min Score (vmin):", value=default_vmin, key="heatmap_vmin", format="%.2f")
#             with col2:
#                 heatmap_vmax = st.number_input("Heatmap Max Score (vmax):", value=default_vmax, key="heatmap_vmax", format="%.2f")

#             if not heatmap_data.empty and not heatmap_data.isna().all().all():
#                 fig_dms, ax_dms = plt.subplots(figsize=(max(12, len(heatmap_data.columns) * 0.7), max(8, len(heatmap_data.index) * 0.7)))
#                 sns.heatmap(
#                     heatmap_data,
#                     cmap='viridis',
#                     annot=False, # Annotate only if heatmap is very small
#                     fmt=".2f",
#                     cbar_kws={'label': 'Average DMS Score'},
#                     vmin=heatmap_vmin,
#                     vmax=heatmap_vmax,
#                     ax=ax_dms
#                 )
#                 ax_dms.set_title(f'Average DMS Score Heatmap (Positions {selected_pos_range[0]} - {selected_pos_range[1]})')
#                 ax_dms.set_xlabel('Protein Position')
#                 ax_dms.set_ylabel('New Amino Acid')
#                 plt.tight_layout()
#                 st.pyplot(fig_dms)
#             else:
#                 st.info("No data available for the selected position range to generate the heatmap. Adjust the range or check data availability in the preview.")

#             # --- Interpretation ---
#             st.subheader("Interpretation: Easiest & Hardest to Mutate Positions")
#             st.write("*(Assuming a **higher** DMS score indicates **easier/more tolerated** mutation, and a **lower** score indicates **harder/less tolerated** mutation)*")

#             # Calculate average score per position (only from filtered data)
#             avg_score_per_position = filtered_dms_df.groupby('position')['DMS_score'].mean().sort_values(ascending=False)
            
#             # Calculate average score per new amino acid (only from filtered data)
#             avg_score_per_new_aa = filtered_dms_df.groupby('new_aa')['DMS_score'].mean().sort_values(ascending=False)

#             col_pos, col_aa = st.columns(2)

#             with col_pos:
#                 if not avg_score_per_position.empty:
#                     st.write("#### Positions by Overall Mutability (Average Score at Position)")
#                     st.write(f"**Easiest to Mutate Positions (Highest Average Score):**")
#                     st.dataframe(avg_score_per_position.head(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                     st.write(f"**Hardest to Mutate Positions (Lowest Average Score):**")
#                     st.dataframe(avg_score_per_position.tail(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                 else:
#                     st.info("No position-wise averages available for the selected range.")

#             with col_aa:
#                 if not avg_score_per_new_aa.empty:
#                     st.write("#### New Amino Acids by General Tolerability (Average Score to New AA)")
#                     st.write(f"**Most Tolerated Mutations (to AA with highest avg score):**")
#                     st.dataframe(avg_score_per_new_aa.head(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                     st.write(f"**Least Tolerated Mutations (to AA with lowest avg score):**")
#                     st.dataframe(avg_score_per_new_aa.tail(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                 else:
#                     st.info("No new amino acid-wise averages available for the selected range.")

#     else:
#         st.info("DMS data is not loaded or is empty. Please check the URL or parsing in the code.")
















# import streamlit as st
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import itertools
# import re # Import regex for parsing mutant string

# # --- Peptide Similarity Calculation Function (re-use from previous code) ---
# def calculate_peptide_similarity(peptide_1: str, peptide_2: str) -> float:
#     peptide_1 = str(peptide_1).strip()
#     peptide_2 = str(peptide_2).strip()

#     if not peptide_1 or not peptide_2:
#         return 0.0

#     if len(peptide_1) != len(peptide_2):
#         return 0.0
            
#     num_matched = 0
#     for i in range(len(peptide_1)):
#         if peptide_1[i] == peptide_2[i]:
#             num_matched += 1
            
#     return num_matched / len(peptide_1)

# # --- Function to calculate similarity scores for a reference peptide (re-use) ---
# @st.cache_data
# def calculate_similarity_scores(reference_peptide: str, df: pd.DataFrame) -> pd.DataFrame:
#     df_copy = df.copy()
#     df_copy['similarity_score'] = float('nan')
    
#     reference_peptide = str(reference_peptide).strip()
#     if not reference_peptide:
#         return df_copy

#     for index, row in df_copy.iterrows():
#         peptide_x = str(row['peptide']).strip()
#         if not peptide_x:
#             continue
#         try:
#             score = calculate_peptide_similarity(reference_peptide, peptide_x)
#             df_copy.at[index, 'similarity_score'] = score
#         except (ValueError, TypeError):
#             continue
#     return df_copy

# # --- Function to determine most similar groups (re-use) ---
# @st.cache_data
# def determine_most_similar_groups(df: pd.DataFrame) -> dict:
#     unique_groups = df['group'].dropna().unique().tolist()
#     if len(unique_groups) < 2:
#         return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}

#     group_pairs = list(itertools.combinations(unique_groups, 2))
#     all_group_pair_scores = {}

#     for group1_name, group2_name in group_pairs:
#         peptides_group1 = df[df['group'] == group1_name]['peptide'].dropna().tolist()
#         peptides_group2 = df[df['group'] == group2_name]['peptide'].dropna().tolist()
#         pairwise_scores_for_current_pair = []

#         if not peptides_group1 or not peptides_group2:
#             all_group_pair_scores[(group1_name, group2_name)] = 0.0
#             continue

#         for p1 in peptides_group1:
#             if not isinstance(p1, str) or not p1.strip(): continue
#             for p2 in peptides_group2:
#                 if not isinstance(p2, str) or not p2.strip(): continue
#                 try:
#                     score = calculate_peptide_similarity(p1.strip(), p2.strip())
#                     pairwise_scores_for_current_pair.append(score)
#                 except (ValueError, TypeError):
#                     continue
        
#         if pairwise_scores_for_current_pair:
#             avg_score = sum(pairwise_scores_for_current_pair) / len(pairwise_scores_for_current_pair)
#             all_group_pair_scores[(group1_name, group2_name)] = avg_score
#         else:
#             all_group_pair_scores[(group1_name, group2_name)] = 0.0

#     if not all_group_pair_scores:
#         return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}
        
#     most_similar_pair = max(all_group_pair_scores.items(), key=lambda item: item[1])

#     return {
#         'group_pair': most_similar_pair[0],
#         'average_similarity_score': most_similar_pair[1],
#         'all_pair_scores': all_group_pair_scores
#     }

# # --- Data Loading Functions ---
# @st.cache_data
# def load_peptide_groups_data():
#     try:
#         df_loaded = pd.read_csv("https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/peptides/peptide_groups.csv")
#         if 'group' not in df_loaded.columns or 'peptide' not in df_loaded.columns:
#             st.error("Error: 'peptide_groups.csv' must contain 'group' and 'peptide' columns.")
#             return pd.DataFrame()
#         df_loaded['peptide'] = df_loaded['peptide'].astype(str).str.strip()
#         df_loaded = df_loaded[df_loaded['peptide'] != ''].copy()
#         df_loaded['group'] = df_loaded['group'].astype(str).str.strip()
#         df_loaded = df_loaded[df_loaded['group'] != ''].copy()
#         return df_loaded
#     except Exception as e:
#         st.error(f"Error loading default 'peptide_groups.csv' from URL: {e}")
#         st.warning("Using a dummy 'Peptide Groups' DataFrame for demonstration due to loading failure.")
#         data = {
#             'group': ['group_A']*50 + ['group_B']*50 + ['group_C']*50,
#             'peptide': [
#                 'AYLDVFNKF', 'VYIDVFNWD', 'AALDVFNKF', 'PYLDVFNKW', 'AYLDFNKKF',
#             ]*10 + [
#                 'VYIDVFNWD', 'GYIDVFNWD', 'TYLDVFNWS', 'VYIDVENWA', 'AYLDVFNKE',
#             ]*10 + [
#                 'PPLDVFNKA', 'QYLDVFNKC', 'RYIDVFNKD', 'AYLDVFNKG', 'AYLDVFNKH'
#             ]*10
#         }
#         import random
#         import string
#         def generate_random_peptide(length=9):
#             return ''.join(random.choice(string.ascii_uppercase) for _ in range(length))
#         for i in range(len(data['peptide'])):
#             if random.random() < 0.2:
#                 data['peptide'][i] = generate_random_peptide()
#         return pd.DataFrame(data)

# # Helper for parsing mutant string
# mutant_pattern = re.compile(r'([A-Z])(\d+)([A-Z])')

# def parse_mutant_column(df, mutant_col='mutant'):
#     parsed_results = df[mutant_col].apply(lambda x: mutant_pattern.match(str(x)))
#     df['original_aa'] = [m.group(1) if m else None for m in parsed_results]
#     df['position'] = [int(m.group(2)) if m and m.group(2).isdigit() else None for m in parsed_results]
#     df['new_aa'] = [m.group(3) if m else None for m in parsed_results]
#     return df

# @st.cache_data
# def load_dms_data():
#     try:
#         dms_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/dms.csv')
#         if 'mutant' not in dms_df.columns or 'DMS_score' not in dms_df.columns:
#             st.error("Error: 'dms.csv' must contain 'mutant' and 'DMS_score' columns.")
#             return pd.DataFrame()
        
#         dms_df = parse_mutant_column(dms_df.copy())
#         dms_df.dropna(subset=['original_aa', 'position', 'new_aa', 'DMS_score'], inplace=True)
#         dms_df['position'] = dms_df['position'].astype(int)
        
#         return dms_df
#     except Exception as e:
#         st.error(f"Error loading 'dms.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
#         return pd.DataFrame()

# @st.cache_data
# def load_if1_data():
#     try:
#         if1_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/if1.csv')
#         if 'mutant' not in if1_df.columns or 'esmif1_ll' not in if1_df.columns:
#             st.error("Error: 'if1.csv' must contain 'mutant' and 'esmif1_ll' columns.")
#             return pd.DataFrame()
        
#         if1_df = parse_mutant_column(if1_df.copy())
#         if1_df.dropna(subset=['original_aa', 'position', 'new_aa', 'esmif1_ll'], inplace=True)
#         if1_df['position'] = if1_df['position'].astype(int)
        
#         return if1_df
#     except Exception as e:
#         st.error(f"Error loading 'if1.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
#         return pd.DataFrame()

# @st.cache_data
# def load_mpnn_data():
#     try:
#         mpnn_df = pd.read_csv('https://raw.githubusercontent.com/Insmed-Computaional-Biology/bi_code_interview/refs/heads/main/proteingym/mpnn.csv')
#         if 'mutant' not in mpnn_df.columns or 'pmpnn_ll' not in mpnn_df.columns:
#             st.error("Error: 'mpnn.csv' must contain 'mutant' and 'pmpnn_ll' columns.")
#             return pd.DataFrame()
        
#         mpnn_df = parse_mutant_column(mpnn_df.copy())
#         mpnn_df.dropna(subset=['original_aa', 'position', 'new_aa', 'pmpnn_ll'], inplace=True)
#         mpnn_df['position'] = mpnn_df['position'].astype(int)
        
#         return mpnn_df
#     except Exception as e:
#         st.error(f"Error loading 'mpnn.csv' from URL: {e}. Please ensure the 'mutant' column format is 'OriginalAA_Position_NewAA' (e.g., A123G).")
#         return pd.DataFrame()


# # --- Streamlit App Layout ---
# st.set_page_config(layout="wide")
# st.title("Protein Analysis Dashboard")
# st.write("Explore peptide group similarities, deep mutational scanning (DMS) data, and compare predictions.")

# # --- Tabbed Interface for different analyses ---
# tab_peptide, tab_dms = st.tabs(["Peptide Group Analysis", "Deep Mutational Scanning (DMS) Analysis"])

# # --- Peptide Group Analysis Tab ---
# with tab_peptide:
#     st.header("Peptide Group Similarity Dashboard")
#     st.write("This dashboard allows you to explore peptide group similarities based on sequence matching.")

#     # --- Data Source Selection for Peptide Analysis ---
#     data_source_option = st.radio(
#         "Choose data source for Peptide Analysis:",
#         ("Use existing 'Peptide Groups' data", "Upload a custom CSV file"),
#         key="peptide_data_source_radio"
#     )

#     df_peptide = None
#     if data_source_option == "Use existing 'Peptide Groups' data":
#         df_peptide = load_peptide_groups_data()
#         if not df_peptide.empty:
#             st.success("Default 'Peptide Groups' data loaded successfully!")
            
#     elif data_source_option == "Upload a custom CSV file":
#         uploaded_file = st.file_uploader("Upload your own CSV file for Peptide Analysis", type="csv", key="peptide_file_uploader")
#         if uploaded_file is not None:
#             try:
#                 df_peptide = pd.read_csv(uploaded_file)
#                 if 'group' not in df_peptide.columns or 'peptide' not in df_peptide.columns:
#                     st.warning("Warning: Uploaded CSV does not contain 'group' and 'peptide' columns. Some features may not work.")
                
#                 df_peptide['peptide'] = df_peptide['peptide'].astype(str).str.strip()
#                 df_peptide = df_peptide[df_peptide['peptide'] != ''].copy()
#                 df_peptide['group'] = df_peptide['group'].astype(str).str.strip()
#                 df_peptide = df_peptide[df_peptide['group'] != ''].copy()

#                 st.success("Custom CSV file uploaded and processed for Peptide Analysis!")

#             except Exception as e:
#                 st.error(f"Error reading uploaded file: {e}. Please ensure it's a valid CSV.")
#                 df_peptide = None
#         else:
#             st.info("Please upload a CSV file to proceed with Peptide Analysis.")


#     # --- Main Peptide Dashboard Logic (Conditional on Data Loading) ---
#     if df_peptide is not None and not df_peptide.empty:
#         st.sidebar.header("Peptide Data Overview")
#         st.sidebar.write(f"**Total Peptides:** {len(df_peptide)}")
#         st.sidebar.write(f"**Number of Groups:** {df_peptide['group'].nunique()}")
        
#         st.subheader("Peptide Data Preview")
#         st.dataframe(df_peptide.head())

#         # --- Section 1: Analyze Reference Peptide Similarity ---
#         st.header("1. Reference Peptide Similarity Analysis")
#         st.write("Determine how similar a chosen peptide is to all peptides within each group.")

#         st.sidebar.subheader("Reference Peptide Configuration")

#         ref_peptide_source_option = st.sidebar.radio(
#             "How to specify Reference Peptide?",
#             ("Select from existing peptides", "Enter custom peptide"),
#             key="ref_peptide_source_radio"
#         )

#         selected_reference_peptide = None 

#         all_peptides_str = "".join(df_peptide['peptide'].dropna().tolist())
#         peptide_alphabet = sorted(list(set(char.upper() for char in all_peptides_str if char.isalpha())))
        
#         if peptide_alphabet:
#             st.sidebar.info(f"**Allowed peptide characters:** {' '.join(peptide_alphabet)}")
#         else:
#             st.sidebar.warning("Could not determine a peptide alphabet from loaded peptide data. Please ensure 'peptide' column contains valid sequences.")


#         if ref_peptide_source_option == "Select from existing peptides":
#             example_peptides = df_peptide['peptide'].dropna().unique().tolist()
#             if not example_peptides:
#                 st.sidebar.warning("No peptides found in data to select from.")
#                 selected_reference_peptide = None
#             else:
#                 default_ref_peptide_val = 'AYLDVFNKF'
#                 default_index = example_peptides.index(default_ref_peptide_val) if default_ref_peptide_val in example_peptides else 0
                
#                 selected_reference_peptide = st.sidebar.selectbox(
#                     "Choose from existing peptides:",
#                     options=example_peptides,
#                     index=default_index,
#                     key="select_existing_ref_peptide"
#                 )

#         elif ref_peptide_source_option == "Enter custom peptide":
#             max_peptide_len = df_peptide['peptide'].str.len().max() if not df_peptide.empty else 20

#             custom_peptide_input = st.sidebar.text_input(
#                 "Enter your custom peptide (e.g., 'AYLDVFNKF'):",
#                 value='AYLDVFNKF', # Default value for custom input
#                 max_chars=max_peptide_len,
#                 key="custom_ref_peptide_input"
#             )
#             custom_peptide_input_clean = custom_peptide_input.strip().upper()

#             if custom_peptide_input_clean:
#                 if not custom_peptide_input_clean.isalpha():
#                     st.sidebar.error("Custom peptide must contain only alphabetical characters (A-Z, a-z).")
#                     selected_reference_peptide = None
#                 else:
#                     is_valid_chars = True
#                     invalid_chars = []
#                     for char in custom_peptide_input_clean:
#                         if char not in peptide_alphabet:
#                             is_valid_chars = False
#                             invalid_chars.append(char)
                    
#                     if not is_valid_chars:
#                         st.sidebar.error(f"Invalid characters detected: {' '.join(set(invalid_chars))}. Only these are allowed: {' '.join(peptide_alphabet)}")
#                         selected_reference_peptide = None
#                     else:
#                         selected_reference_peptide = custom_peptide_input_clean
#             else:
#                 st.sidebar.info("Please enter a peptide sequence.")
#                 selected_reference_peptide = None

#         all_groups = df_peptide['group'].dropna().unique().tolist()
#         if not all_groups:
#             st.warning("No peptide groups found in the loaded data to analyze.")
#         elif selected_reference_peptide:
#             selected_groups_for_ref = st.multiselect(
#                 "Select Groups to Analyze (for reference peptide):",
#                 options=all_groups,
#                 default=all_groups,
#                 key="ref_peptide_group_multiselect"
#             )

#             if selected_groups_for_ref:
#                 df_filtered_for_ref = df_peptide[df_peptide['group'].isin(selected_groups_for_ref)]
                
#                 df_with_scores = calculate_similarity_scores(selected_reference_peptide, df_filtered_for_ref)
#                 valid_scores_df = df_with_scores.dropna(subset=['similarity_score'])

#                 if not valid_scores_df.empty:
#                     st.subheader(f"Similarity Score Distribution to '{selected_reference_peptide}'")
#                     fig1, ax1 = plt.subplots(figsize=(12, 7))
#                     sns.boxplot(x='group', y='similarity_score', data=valid_scores_df, ax=ax1, palette='viridis')
#                     ax1.set_title(f'Similarity Score Distribution to "{selected_reference_peptide}" Across Selected Groups')
#                     ax1.set_xlabel('Peptide Group')
#                     ax1.set_ylabel(f'Similarity Score to "{selected_reference_peptide}"')
#                     ax1.set_ylim(0, 1)
#                     ax1.grid(axis='y', linestyle='--', alpha=0.7)
#                     st.pyplot(fig1)

#                     st.subheader("Quantitative Spread within Each Group:")
#                     group_spread_stats = valid_scores_df.groupby('group')['similarity_score'].agg(
#                         mean_score='mean',
#                         median_score='median',
#                         std_dev_score='std',
#                         iqr_score=lambda x: x.quantile(0.75) - x.quantile(0.25),
#                         count='count'
#                     ).dropna()
#                     st.dataframe(group_spread_stats.style.format("{:.4f}"))
#                 else:
#                     st.info(f"No valid similarity scores found for '{selected_reference_peptide}' within the selected groups. This might be due to peptide length mismatch or group content.")
#             else:
#                 st.info("Please select at least one group to analyze against the reference peptide.")
#         else:
#             st.info("Please select or enter a valid reference peptide to proceed with analysis.")


#         # --- Section 2: Explore Inter-Group Similarity ---
#         st.header("2. Inter-Group Similarity Analysis")
#         st.write("Identify and visualize how similar different peptide groups are to each other.")

#         group_similarity_results = determine_most_similar_groups(df_peptide.copy())

#         if group_similarity_results['group_pair'] and group_similarity_results['all_pair_scores']:
#             sorted_pair_scores_list = sorted(group_similarity_results['all_pair_scores'].items(), key=lambda item: item[1], reverse=True)
#             pair_scores_df = pd.DataFrame(
#                 [(f"{pair[0]} vs {pair[1]}", score) for pair, score in sorted_pair_scores_list],
#                 columns=['Group Pair', 'Average Similarity Score']
#             )

#             st.subheader("Average Similarity Scores Between All Group Pairs:")
#             num_top_pairs = st.slider("Show Top N Most Similar Pairs:", 1, len(pair_scores_df), min(5, len(pair_scores_df)), key="num_top_pairs_slider")
            
#             fig2, ax2 = plt.subplots(figsize=(12, max(6, num_top_pairs * 0.8)))
#             sns.barplot(x='Average Similarity Score', y='Group Pair', data=pair_scores_df.head(num_top_pairs), palette='coolwarm', ax=ax2)
#             ax2.set_title(f'Top {num_top_pairs} Most Similar Peptide Group Pairs')
#             ax2.set_xlabel('Average Similarity Score (0-1)')
#             ax2.set_ylabel('Peptide Group Pair')
#             ax2.set_xlim(0, 1)
#             ax2.grid(axis='x', linestyle='--', alpha=0.7)
#             st.pyplot(fig2)

#             st.success(f"**The single most similar pair is: {group_similarity_results['group_pair'][0]} and {group_similarity_results['group_pair'][1]}** (Average Score: {group_similarity_results['average_similarity_score']:.4f})")

#             st.subheader("Detailed Pairwise Similarity Heatmap")
            
#             available_pairs_for_heatmap = [f"{pair[0]} vs {pair[1]}" for pair, score in sorted_pair_scores_list]
#             selected_heatmap_pair_str = st.selectbox(
#                 "Select a Group Pair for Detailed Heatmap:",
#                 options=available_pairs_for_heatmap,
#                 index=0,
#                 key="heatmap_pair_selectbox"
#             )

#             selected_group1_name, selected_group2_name = selected_heatmap_pair_str.split(' vs ')

#             peptides_g1 = df_peptide[df_peptide['group'] == selected_group1_name]['peptide'].dropna().tolist()
#             peptides_g2 = df_peptide[df_peptide['group'] == selected_group2_name]['peptide'].dropna().tolist()

#             max_heatmap_peptides = st.slider(
#                 "Limit # of peptides for heatmap (adjust for clarity):",
#                 min_value=5, max_value=min(50, max(len(peptides_g1), len(peptides_g2), 10)),
#                 value=min(20, max(len(peptides_g1), len(peptides_g2), 10)),
#                 key="heatmap_peptide_limit"
#             )
#             peptides_g1_sampled = peptides_g1[:max_heatmap_peptides]
#             peptides_g2_sampled = peptides_g2[:max_heatmap_peptides]

#             similarity_matrix_inter_group = pd.DataFrame(
#                 np.nan,
#                 index=peptides_g1_sampled,
#                 columns=peptides_g2_sampled,
#                 dtype=float
#             )

#             for p1 in peptides_g1_sampled:
#                 for p2 in peptides_g2_sampled:
#                     if not isinstance(p1, str) or not p1.strip() or not isinstance(p2, str) or not p2.strip():
#                         similarity_matrix_inter_group.loc[p1, p2] = np.nan
#                         continue
#                     try:
#                         score = calculate_peptide_similarity(p1, p2)
#                         similarity_matrix_inter_group.loc[p1, p2] = score
#                     except (ValueError, TypeError):
#                         similarity_matrix_inter_group.loc[p1, p2] = np.nan
#                         continue

#             if not similarity_matrix_inter_group.empty and not similarity_matrix_inter_group.isna().all().all():
#                 fig4, ax4 = plt.subplots(figsize=(max(10, len(peptides_g2_sampled)*0.8), max(8, len(peptides_g1_sampled)*0.8)))
#                 sns.heatmap(
#                     similarity_matrix_inter_group,
#                     cmap='viridis',
#                     annot=True if (len(peptides_g1_sampled) < 15 and len(peptides_g2_sampled) < 15) else False,
#                     fmt=".2f",
#                     cbar_kws={'label': 'Similarity Score (0-1)'},
#                     ax=ax4
#                 )
#                 ax4.set_title(f'Pairwise Similarity Between {selected_group1_name} and {selected_group2_name}')
#                 ax4.set_xlabel(f'Peptides in {selected_group2_name}')
#                 ax4.set_ylabel(f'Peptides in {selected_group1_name}')
#                 plt.tight_layout()
#                 st.pyplot(fig4)
#             else:
#                 st.info("No valid data or peptides for this heatmap. Adjust selections or check data content (e.g., peptide strings might be empty).")
#         else:
#             st.warning("Not enough distinct peptide groups to perform inter-group similarity analysis (need at least 2), or no valid peptide data loaded.")

#     else:
#         st.info("Please load peptide data using the options above to enable analysis.")


# # --- Deep Mutational Scanning (DMS) Analysis Tab (now includes Prediction Comparison) ---
# with tab_dms:
#     st.header("Deep Mutational Scanning (DMS) Analysis & Prediction Comparison")
#     st.write("Visualize DMS scores as a heatmap to identify positions and mutations affecting protein fitness/function, and compare them with predicted scores.")

#     dms_df_main = load_dms_data() # Main DMS data for this tab
#     if1_df_main = load_if1_data() # Main IF1 data for this tab
#     mpnn_df_main = load_mpnn_data() # New: Main MPNN data for this tab

#     # --- Section 1: DMS Heatmap and Interpretation ---
#     st.subheader("1. DMS Score Heatmap Analysis")
#     if not dms_df_main.empty:
#         st.success("DMS data loaded successfully!")
#         st.subheader("DMS Data Preview")
#         st.dataframe(dms_df_main.head())

#         st.subheader("DMS Heatmap Configuration")

#         all_positions = sorted(dms_df_main['position'].unique())
#         all_new_aas = sorted(dms_df_main['new_aa'].unique())

#         if not all_positions or not all_new_aas:
#             st.warning("Not enough data to generate heatmap (missing positions or new amino acids after parsing).")
#         else:
#             min_pos, max_pos = min(all_positions), max(all_positions)
#             selected_pos_range = st.slider(
#                 "Select Position Range for Heatmap:",
#                 min_value=min_pos,
#                 max_value=max_pos,
#                 value=(min_pos, max_pos),
#                 key="dms_pos_range_slider"
#             )

#             filtered_dms_df = dms_df_main[(dms_df_main['position'] >= selected_pos_range[0]) & 
#                                      (dms_df_main['position'] <= selected_pos_range[1])]

#             heatmap_data = filtered_dms_df.pivot_table(
#                 values='DMS_score',
#                 index='new_aa',
#                 columns='position',
#                 aggfunc='mean'
#             ).fillna(np.nan)

#             heatmap_data = heatmap_data.reindex(index=all_new_aas, columns=range(min_pos, max_pos + 1), fill_value=np.nan)
            
#             st.subheader("DMS Score Heatmap: Average Score per Mutation")
#             st.write("Colors represent the average DMS score for mutations. Brighter (yellow) indicates higher scores, suggesting mutations are more tolerated at that position to that new amino acid.")

#             valid_scores_flat = heatmap_data.stack().dropna()
#             if not valid_scores_flat.empty:
#                 default_vmin = valid_scores_flat.min()
#                 default_vmax = valid_scores_flat.max()
#             else:
#                 default_vmin = 0.0
#                 default_vmax = 1.0

#             col1, col2 = st.columns(2)
#             with col1:
#                 heatmap_vmin = st.number_input("Heatmap Min Score (vmin):", value=default_vmin, key="heatmap_vmin", format="%.2f")
#             with col2:
#                 heatmap_vmax = st.number_input("Heatmap Max Score (vmax):", value=default_vmax, key="heatmap_vmax", format="%.2f")

#             if not heatmap_data.empty and not heatmap_data.isna().all().all():
#                 fig_dms, ax_dms = plt.subplots(figsize=(max(12, len(heatmap_data.columns) * 0.7), max(8, len(heatmap_data.index) * 0.7)))
#                 sns.heatmap(
#                     heatmap_data,
#                     cmap='viridis',
#                     annot=False,
#                     fmt=".2f",
#                     cbar_kws={'label': 'Average DMS Score'},
#                     vmin=heatmap_vmin,
#                     vmax=heatmap_vmax,
#                     ax=ax_dms
#                 )
#                 ax_dms.set_title(f'Average DMS Score Heatmap (Positions {selected_pos_range[0]} - {selected_pos_range[1]})')
#                 ax_dms.set_xlabel('Protein Position')
#                 ax_dms.set_ylabel('New Amino Acid')
#                 plt.tight_layout()
#                 st.pyplot(fig_dms)
#             else:
#                 st.info("No data available for the selected position range to generate the heatmap. Adjust the range or check data availability in the preview.")

#             st.subheader("Interpretation: Easiest & Hardest to Mutate Positions")
#             st.write("*(Assuming a **higher** DMS score indicates **easier/more tolerated** mutation, and a **lower** score indicates **harder/less tolerated** mutation)*")

#             avg_score_per_position = filtered_dms_df.groupby('position')['DMS_score'].mean().sort_values(ascending=False)
#             avg_score_per_new_aa = filtered_dms_df.groupby('new_aa')['DMS_score'].mean().sort_values(ascending=False)

#             col_pos, col_aa = st.columns(2)

#             with col_pos:
#                 if not avg_score_per_position.empty:
#                     st.write("#### Positions by Overall Mutability (Average Score at Position)")
#                     st.write(f"**Easiest to Mutate Positions (Highest Average Score):**")
#                     st.dataframe(avg_score_per_position.head(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                     st.write(f"**Hardest to Mutate Positions (Lowest Average Score):**")
#                     st.dataframe(avg_score_per_position.tail(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                 else:
#                     st.info("No position-wise averages available for the selected range.")

#             with col_aa:
#                 if not avg_score_per_new_aa.empty:
#                     st.write("#### New Amino Acids by General Tolerability (Average Score to New AA)")
#                     st.write(f"**Most Tolerated Mutations (to AA with highest avg score):**")
#                     st.dataframe(avg_score_per_new_aa.head(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                     st.write(f"**Least Tolerated Mutations (to AA with lowest avg score):**")
#                     st.dataframe(avg_score_per_new_aa.tail(5).reset_index().rename(columns={'DMS_score': 'Avg DMS Score'}))
#                 else:
#                     st.info("No new amino acid-wise averages available for the selected range.")
#     else:
#         st.info("DMS data is not loaded or is empty for heatmap analysis. Please check the URL or parsing in the code.")


#     # --- Section 2: Prediction Comparison (Moved from its own tab and expanded) ---
#     st.markdown("---") # Add a separator for visual clarity
#     st.header("2. DMS vs. Prediction Comparison")
#     st.write("Compare experimental Deep Mutational Scanning (DMS) scores with predicted scores from different models (esmif1_ll, pmpnn_ll).")

#     if not dms_df_main.empty and not if1_df_main.empty and not mpnn_df_main.empty:
#         st.success("All DMS, IF1, and MPNN prediction data loaded successfully for comparison!")

#         # Merge for IF1 comparison
#         comparison_df_if1 = pd.merge(dms_df_main[['mutant', 'DMS_score']], 
#                                      if1_df_main[['mutant', 'esmif1_ll']], 
#                                      on='mutant', 
#                                      how='inner')
        
#         # Merge for MPNN comparison
#         comparison_df_mpnn = pd.merge(dms_df_main[['mutant', 'DMS_score']], 
#                                       mpnn_df_main[['mutant', 'pmpnn_ll']], 
#                                       on='mutant', 
#                                       how='inner')
        
#         if comparison_df_if1.empty and comparison_df_mpnn.empty:
#             st.warning("No common mutations found between DMS and either prediction dataset. Please check 'mutant' column formatting and data content.")
#         else:
#             st.subheader("Overall Model Similarity to DMS Data")

#             correlation_if1 = 0.0
#             if not comparison_df_if1.empty:
#                 correlation_if1 = comparison_df_if1['DMS_score'].corr(comparison_df_if1['esmif1_ll'])

#             correlation_mpnn = 0.0
#             if not comparison_df_mpnn.empty:
#                 correlation_mpnn = comparison_df_mpnn['DMS_score'].corr(comparison_df_mpnn['pmpnn_ll'])

#             col_corr1, col_corr2 = st.columns(2)
#             with col_corr1:
#                 st.metric(label="Pearson Correlation (DMS vs. IF1)", value=f"{correlation_if1:.4f}")
#                 if comparison_df_if1.empty:
#                     st.info("No common mutations for IF1 correlation.")
#             with col_corr2:
#                 st.metric(label="Pearson Correlation (DMS vs. MPNN)", value=f"{correlation_mpnn:.4f}")
#                 if comparison_df_mpnn.empty:
#                     st.info("No common mutations for MPNN correlation.")
            
#             st.markdown("---") # Visual separator

#             if not comparison_df_if1.empty or not comparison_df_mpnn.empty:
#                 if correlation_if1 > correlation_mpnn:
#                     st.success(f"**The IF1 model (Pearson R: {correlation_if1:.4f}) is more similar to DMS data than the MPNN model (Pearson R: {correlation_mpnn:.4f}).**")
#                 elif correlation_mpnn > correlation_if1:
#                     st.success(f"**The MPNN model (Pearson R: {correlation_mpnn:.4f}) is more similar to DMS data than the IF1 model (Pearson R: {correlation_if1:.4f}).**")
#                 else:
#                     st.info("The IF1 and MPNN models have very similar correlation values with DMS data.")
#             else:
#                 st.warning("Not enough data to determine which model is more similar (no common mutations with DMS).")


#             st.markdown("---") # Another visual separator
#             st.subheader("Detailed Comparison for Selected Model")

#             model_to_view = st.selectbox(
#                 "Select Prediction Model to View Details:",
#                 options=['IF1 (esmif1_ll)', 'MPNN (pmpnn_ll)'],
#                 key="model_comparison_selector"
#             )

#             if model_to_view == 'IF1 (esmif1_ll)':
#                 current_comparison_df = comparison_df_if1
#                 predicted_score_col = 'esmif1_ll'
#                 model_name = 'IF1'
#             else: # MPNN (pmpnn_ll)
#                 current_comparison_df = comparison_df_mpnn
#                 predicted_score_col = 'pmpnn_ll'
#                 model_name = 'MPNN'

#             if not current_comparison_df.empty:
#                 # Scatter Plot for selected model
#                 st.write(f"#### Scatter Plot: Experimental vs. Predicted Scores for {model_name}")
#                 st.write("Each point represents a mutation. Closeness to the diagonal line (y=x) indicates better agreement.")

#                 fig_scatter, ax_scatter = plt.subplots(figsize=(10, 8))
#                 sns.scatterplot(x='DMS_score', y=predicted_score_col, data=current_comparison_df, ax=ax_scatter, alpha=0.6)
                
#                 max_val = max(current_comparison_df['DMS_score'].max(), current_comparison_df[predicted_score_col].max())
#                 min_val = min(current_comparison_df['DMS_score'].min(), current_comparison_df[predicted_score_col].min())
#                 ax_scatter.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', label='Perfect Agreement')

#                 ax_scatter.set_title(f'DMS Score vs. {predicted_score_col} ({model_name} Predicted Score)')
#                 ax_scatter.set_xlabel('DMS Score (Experimental)')
#                 ax_scatter.set_ylabel(f'{predicted_score_col} (Predicted)')
#                 ax_scatter.legend()
#                 ax_scatter.grid(True, linestyle='--', alpha=0.7)
#                 st.pyplot(fig_scatter)

#                 current_comparison_df['absolute_difference'] = abs(current_comparison_df['DMS_score'] - current_comparison_df[predicted_score_col])

#                 # Highest Agreement (Lowest Difference) for selected model
#                 st.write(f"#### Mutations with Highest Agreement ({model_name} Model)")
#                 st.write("These mutations have the smallest absolute difference between experimental and predicted scores.")
#                 highest_agreement_mutations = current_comparison_df.sort_values(by='absolute_difference', ascending=True).head(10)
#                 st.dataframe(highest_agreement_mutations[['mutant', 'DMS_score', predicted_score_col, 'absolute_difference']].style.format({"DMS_score": "{:.4f}", predicted_score_col: "{:.4f}", "absolute_difference": "{:.4f}"}))

#                 # Lowest Agreement (Highest Difference) for selected model
#                 st.write(f"#### Mutations with Lowest Agreement ({model_name} Model)")
#                 st.write("These mutations have the largest absolute difference between experimental and predicted scores, indicating significant disagreement.")
#                 lowest_agreement_mutations = current_comparison_df.sort_values(by='absolute_difference', ascending=False).head(10)
#                 st.dataframe(lowest_agreement_mutations[['mutant', 'DMS_score', predicted_score_col, 'absolute_difference']].style.format({"DMS_score": "{:.4f}", predicted_score_col: "{:.4f}", "absolute_difference": "{:.4f}"}))
#             else:
#                 st.info(f"No common mutations found for {model_name} to perform detailed comparison.")

#     elif dms_df_main.empty:
#         st.warning("DMS data not available for prediction comparison. Please ensure DMS data loads successfully.")
#     elif if1_df_main.empty:
#         st.warning("Prediction data (if1.csv) not loaded. Cannot perform prediction comparison.")
#     elif mpnn_df_main.empty:
#         st.warning("Prediction data (mpnn.csv) not loaded. Cannot perform prediction comparison.")
#     else:
#         st.info("Loading all data for prediction comparison...")





import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import re # Import regex for parsing mutant string

# --- Peptide Similarity Calculation Function (re-use from previous code) ---
def calculate_peptide_similarity(peptide_1: str, peptide_2: str) -> float:
    peptide_1 = str(peptide_1).strip()
    peptide_2 = str(peptide_2).strip()

    if not peptide_1 or not peptide_2:
        return 0.0

    if len(peptide_1) != len(peptide_2):
        return 0.0
            
    num_matched = 0
    for i in range(len(peptide_1)):
        if peptide_1[i] == peptide_2[i]:
            num_matched += 1
            
    return num_matched / len(peptide_1)

# --- Function to calculate similarity scores for a reference peptide (re-use) ---
@st.cache_data
def calculate_similarity_scores(reference_peptide: str, df: pd.DataFrame) -> pd.DataFrame:
    df_copy = df.copy()
    df_copy['similarity_score'] = float('nan')
    
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

# --- Function to determine most similar groups (re-use) ---
@st.cache_data
def determine_most_similar_groups(df: pd.DataFrame) -> dict:
    unique_groups = df['group'].dropna().unique().tolist()
    if len(unique_groups) < 2:
        return {'group_pair': None, 'average_similarity_score': -1.0, 'all_pair_scores': {}}

    group_pairs = list(itertools.combinations(unique_groups, 2))
    all_group_pair_scores = {}

    for group1_name, group2_name in group_pairs:
        peptides_group1 = df[df['group'] == group1_name]['peptide'].dropna().tolist()
        peptides_group2 = df[df['group'] == group2_name]['peptide'].dropna().tolist()
        pairwise_scores_for_current_pair = []

        if not peptides_group1 or not peptides_group2:
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
        
    most_similar_pair = max(all_group_pair_scores.items(), key=lambda item: item[1])

    return {
        'group_pair': most_similar_pair[0],
        'average_similarity_score': most_similar_pair[1],
        'all_pair_scores': all_group_pair_scores
    }

# --- Data Loading Functions ---
@st.cache_data
def load_peptide_groups_data():
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
mutant_pattern = re.compile(r'([A-Z])(\d+)([A-Z])')

def parse_mutant_column(df, mutant_col='mutant'):
    parsed_results = df[mutant_col].apply(lambda x: mutant_pattern.match(str(x)))
    df['original_aa'] = [m.group(1) if m else None for m in parsed_results]
    df['position'] = [int(m.group(2)) if m and m.group(2).isdigit() else None for m in parsed_results]
    df['new_aa'] = [m.group(3) if m else None for m in parsed_results]
    return df

@st.cache_data
def load_dms_data():
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


# --- Streamlit App Layout ---
st.set_page_config(layout="wide")
st.title("Protein Analysis Dashboard")
st.write("Explore peptide group similarities, deep mutational scanning (DMS) data, and compare predictions.")

# --- Tabbed Interface for different analyses ---
tab_peptide, tab_dms = st.tabs(["Peptide Group Analysis", "Deep Mutational Scanning (DMS) Analysis"])

# --- Peptide Group Analysis Tab ---
with tab_peptide:
    st.header("Peptide Group Similarity Dashboard")
    st.write("This dashboard allows you to explore peptide group similarities based on sequence matching.")

    # --- Data Source Selection for Peptide Analysis ---
    # Removed the "Upload a custom CSV file" option
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
            
    # Removed the 'elif' block for custom CSV upload

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


        all_groups = df_peptide['group'].dropna().unique().tolist()
        if not all_groups:
            st.warning("No peptide groups found in the loaded data to analyze.")
        elif selected_reference_peptide: # This checks if a reference peptide was successfully determined
            selected_groups_for_ref = st.multiselect(
                "Select Groups to Analyze (for reference peptide):",
                options=all_groups,
                default=all_groups,
                key="ref_peptide_group_multiselect"
            )

            if selected_groups_for_ref:
                df_filtered_for_ref = df_peptide[df_peptide['group'].isin(selected_groups_for_ref)]
                
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
                    ax1.grid(axis='y', linestyle='--', alpha=0.7)
                    st.pyplot(fig1)

                    st.subheader("Quantitative Spread within Each Group:")
                    group_spread_stats = valid_scores_df.groupby('group')['similarity_score'].agg(
                        mean_score='mean',
                        median_score='median',
                        std_dev_score='std',
                        iqr_score=lambda x: x.quantile(0.75) - x.quantile(0.25),
                        count='count'
                    ).dropna()
                    st.dataframe(group_spread_stats.style.format("{:.4f}"))
                else:
                    st.info(f"No valid similarity scores found for '{selected_reference_peptide}' within the selected groups. This might be due to peptide length mismatch or group content.")
            else:
                st.info("Please select at least one group to analyze against the reference peptide.")
        else:
            st.info("Please enter a valid reference peptide (9 characters, allowed characters) to proceed with this analysis section.")


        # --- Section 2: Explore Inter-Group Similarity ---
        # This entire section now only displays if the reference peptide is valid
        if is_reference_peptide_valid:
            st.header("2. Inter-Group Similarity Analysis")
            st.write("Identify and visualize how similar different peptide groups are to each other.")

            group_similarity_results = determine_most_similar_groups(df_peptide.copy())

            if group_similarity_results['group_pair'] and group_similarity_results['all_pair_scores']:
                sorted_pair_scores_list = sorted(group_similarity_results['all_pair_scores'].items(), key=lambda item: item[1], reverse=True)
                pair_scores_df = pd.DataFrame(
                    [(f"{pair[0]} vs {pair[1]}", score) for pair, score in sorted_pair_scores_list],
                    columns=['Group Pair', 'Average Similarity Score']
                )

                st.subheader("Average Similarity Scores Between All Group Pairs:")
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

                max_heatmap_peptides = st.slider(
                    "Limit # of peptides for heatmap (adjust for clarity):",
                    min_value=5, max_value=min(50, max(len(peptides_g1), len(peptides_g2), 10)),
                    value=min(20, max(len(peptides_g1), len(peptides_g2), 10)),
                    key="heatmap_peptide_limit"
                )
                peptides_g1_sampled = peptides_g1[:max_heatmap_peptides]
                peptides_g2_sampled = peptides_g2[:max_heatmap_peptides]

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
                    ax4.set_title(f'Pairwise Similarity Between {selected_group1_name} and {selected_group2_name}')
                    ax4.set_xlabel(f'Peptides in {selected_group2_name}')
                    ax4.set_ylabel(f'Peptides in {selected_group1_name}')
                    plt.tight_layout()
                    st.pyplot(fig4)
                else:
                    st.info("No valid data or peptides for this heatmap. Adjust selections or check data content (e.g., peptide strings might be empty).")
            else:
                st.warning("Not enough distinct peptide groups to perform inter-group similarity analysis (need at least 2), or no valid peptide data loaded.")
        else:
            st.info("Please enter a valid reference peptide (9 characters long, with allowed characters) to enable 'Inter-Group Similarity Analysis'.")

    else:
        st.info("Please load peptide data using the options above to enable analysis.")


# --- Deep Mutational Scanning (DMS) Analysis Tab (now includes Prediction Comparison) ---
with tab_dms:
    st.header("Deep Mutational Scanning (DMS) Analysis & Prediction Comparison")
    st.write("Visualize DMS scores as a heatmap to identify positions and mutations affecting protein fitness/function, and compare them with predicted scores.")

    dms_df_main = load_dms_data() # Main DMS data for this tab
    if1_df_main = load_if1_data() # Main IF1 data for this tab
    mpnn_df_main = load_mpnn_data() # New: Main MPNN data for this tab

    # --- Section 1: DMS Heatmap and Interpretation ---
    st.subheader("1. DMS Score Heatmap Analysis")
    if not dms_df_main.empty:
        st.success("DMS data loaded successfully!")
        st.subheader("DMS Data Preview")
        st.dataframe(dms_df_main.head())

        st.subheader("DMS Heatmap Configuration")

        all_positions = sorted(dms_df_main['position'].unique())
        all_new_aas = sorted(dms_df_main['new_aa'].unique())

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

            filtered_dms_df = dms_df_main[(dms_df_main['position'] >= selected_pos_range[0]) & 
                                     (dms_df_main['position'] <= selected_pos_range[1])]

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

            valid_scores_flat = heatmap_data.stack().dropna()
            if not valid_scores_flat.empty:
                default_vmin = valid_scores_flat.min()
                default_vmax = valid_scores_flat.max()
            else:
                default_vmin = 0.0
                default_vmax = 1.0

            col1, col2 = st.columns(2)
            with col1:
                heatmap_vmin = st.number_input("Heatmap Min Score (vmin):", value=default_vmin, key="heatmap_vmin", format="%.2f")
            with col2:
                heatmap_vmax = st.number_input("Heatmap Max Score (vmax):", value=default_vmax, key="heatmap_vmax", format="%.2f")

            if not heatmap_data.empty and not heatmap_data.isna().all().all():
                # --- Improved Heatmap Plotting ---
                num_aas = len(heatmap_data.index)
                num_positions = len(heatmap_data.columns)

                # Dynamically adjust figure size for better visibility
                fig_dms, ax_dms = plt.subplots(figsize=(max(15, num_positions * 0.7), max(10, num_aas * 0.5))) 
                
                # Conditional annotation based on number of cells
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
            st.markdown(
                """
                <p>
                (Assuming a higher DMS score indicates <span style="border-bottom: 3px solid yellow;">easier/more tolerated</span> mutation, 
                and a lower score indicates <span style="border-bottom: 3px solid purple;">harder/less tolerated</span> mutation)
                </p>
                """, 
                unsafe_allow_html=True
            )

            avg_score_per_position = filtered_dms_df.groupby('position')['DMS_score'].mean().sort_values(ascending=False)
            avg_score_per_new_aa = filtered_dms_df.groupby('new_aa')['DMS_score'].mean().sort_values(ascending=False)

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


    # --- Section 2: Prediction Comparison (Moved from its own tab and expanded) ---
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
        
        if comparison_df_if1.empty and comparison_df_mpnn.empty:
            st.warning("No common mutations found between DMS and either prediction dataset. Please check 'mutant' column formatting and data content.")
        else:
            st.subheader("Overall Model Similarity to DMS Data")

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


            st.markdown("---") # Another visual separator
            st.subheader("Detailed Comparison for Selected Model")

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
                sns.scatterplot(x='DMS_score', y=predicted_score_col, data=current_comparison_df, ax=ax_scatter, alpha=0.6)
                
                max_val = max(current_comparison_df['DMS_score'].max(), current_comparison_df[predicted_score_col].max())
                min_val = min(current_comparison_df['DMS_score'].min(), current_comparison_df[predicted_score_col].min())
                ax_scatter.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', label='Perfect Agreement')

                ax_scatter.set_title(f'DMS Score vs. {predicted_score_col} ({model_name} Predicted Score)')
                ax_scatter.set_xlabel('DMS Score (Experimental)')
                ax_scatter.set_ylabel(f'{predicted_score_col} (Predicted)')
                ax_scatter.legend()
                ax_scatter.grid(True, linestyle='--', alpha=0.7)
                st.pyplot(fig_scatter)

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