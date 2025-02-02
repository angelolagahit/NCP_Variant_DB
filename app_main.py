import streamlit as st
import os
import json
import pandas as pd
import yaml
from variant.IndexHPO import IndexHPO
from variant.HPOAnnotator import HPOAnnotator
from variant.dataframe import *
from variant.query import *
from variant.validator import *



# Paths for storing the HPO index
OUTPUT_FOLDER = "output"
INDEX_FILE = os.path.join(OUTPUT_FOLDER, "hp.index")

# Streamlit UI
st.set_page_config(layout="wide")
st.title("🧬 PhenoVariant App")

# Navigation
tab1, tab2, tab3 = st.tabs(["🧬 Variant Database", "🔍 HPO Annotation", "📚 About"])



# Tab 1: Variant Database
with tab1:
    # Sidebar: Upload YAML Configuration
    st.sidebar.header("🛠 Configuration")
    st.sidebar.subheader("Variant Database")
    yaml_file = st.sidebar.file_uploader("Upload `config.yaml` file", type=["yaml", "yml"])
    config = None

    if yaml_file:
        try:
            config = yaml.safe_load(yaml_file)
            base_path = config['base_path']
            files = config['files']
            st.sidebar.success("✅ YAML configuration loaded successfully!")
        except Exception as e:
            st.sidebar.error(f"❌ Failed to load YAML: {e}")
            config = None

    # Navigation: Variant Database
    options = ["Load Data", "Standardise Data", "Query Data", "Validate Variants"]
    variant_tab = st.radio("Select option:", options, horizontal=False)

    # Initialising standard column names
    if config:
        standard_columns = [
            'MRN', 'Patient Name', 'Var Count', 'Var Number', 'Phenotype', 'Solved Status',
            'Gene', 'HGVSg', 'Transcript', 'Variant (HGVSc)', 'HGVSp', 'Type', 'Zygosity', 'Inheritance',
        ]

        # Stores renamed column names from each dataframe
        renaming_maps = {}

    # Load Data Tab
    if variant_tab == "Load Data":
        st.header("📂 Load Variant Data")

        if config:
            try:
                with st.spinner("Loading data..."):
                    # Reading files
                    lab_cases_df = pd.read_excel(f"{base_path}{files['lab_cases']}", sheet_name='Sheet1', header=0)
                    invitae_summary_df = pd.read_excel(f"{base_path}{files['invitae_summary']}", sheet_name='Invitae list header', header=0)
                    clinical_summary_df = pd.read_excel(f"{base_path}{files['clinical_summary']}", sheet_name='11 Dec', header=0)
                    research_summary_df = pd.read_excel(f"{base_path}{files['research_summary']}", sheet_name='Overall List', header=2)

                    # Handling unconventional ATM summary format
                    atm_summary_df = pd.read_excel(f"{base_path}{files['atm_summary']}", sheet_name='SUMMARY', header=None)
                    atm_summary_df_filtered = atm_summary_df.iloc[8:17].reset_index(drop=True)
                    atm_structured_df = atm_summary_df_filtered.set_index(0).T.reset_index(drop=True)

                    # Store imported dataframes into a dictionary
                    dataframes = {
                        "Lab Cases": lab_cases_df,
                        "ATM Summary": atm_structured_df,
                        "Invitae Summary": invitae_summary_df,
                        "Clinical Summary": clinical_summary_df,
                        "Research Summary": research_summary_df,
                    }

                    st.session_state["dataframes"] = dataframes  # Store in session state
                    st.success("✅ Data successfully loaded!")

                    # Display loaded dataframes
                    for name, df in dataframes.items():
                        st.write(f"### {name}", df.head())

            except FileNotFoundError as e:
                st.error(f"❌ File not found: {e}")
            except ValueError as e:
                st.error(f"❌ Error reading sheet or invalid data format: {e}")
        else:
            st.warning("Upload a YAML config file.")



   # Standardise Data Tab
    elif variant_tab == "Standardise Data":
        st.header("🛠 Standardise Data")

        if st.session_state["dataframes"]:
            try:
                with st.spinner("Standardising and combining data..."):
                    dataframes = st.session_state["dataframes"]
                    
                    # Standardising Lab Cases
                    lab_cases_cols = {
                        'G4K Sample ID': 'MRN',
                        'Number variants detected': 'Var Count',
                        'Variant_Number': 'Var Number',
                        'gene': 'Gene',
                        'type': 'Type',
                        'zygosity': 'Zygosity',
                        'inheritance': 'Inheritance'
                    }

                    std_lab_cases_df = standardise_and_transform_lab_cases(dataframes["Lab Cases"])
                    std_lab_cases_df.rename(columns=lab_cases_cols, inplace=True)
                    std_lab_cases_df = std_lab_cases_df.loc[:, std_lab_cases_df.columns.intersection(standard_columns)]
                    std_lab_cases_df = std_lab_cases_df.reindex(columns=standard_columns)
                    remove_duplicate_phenotypes(std_lab_cases_df)
                    std_lab_cases_df = filter_invalid_variants(std_lab_cases_df)

                    
                    # Standardising Invitae Summary
                    invitae_summary_cols = {
                        'Patient ID (MRN)': 'MRN',
                        'Patient Name': 'Patient Name',
                        'Result': 'Solved Status',
                        'Gene': 'Gene',
                        'Transcript': 'Transcript',
                        'Variant': 'Variant (HGVSc)',
                        'Protein Change': 'HGVSp',
                    }
                    
                    std_invitae_summary_df = standardise_dataframe(dataframes["Invitae Summary"], invitae_summary_cols, standard_columns)
                    
                    # Standardising Clinical Summary
                    clinical_summary_cols = {
                        'Identification No.': 'MRN',
                        'Medical Prob description': 'Phenotype',
                        'Patient Name': 'Patient Name'
                    }
                    
                    std_clinical_summary_df = standardise_dataframe(dataframes["Clinical Summary"], clinical_summary_cols, standard_columns)
                    
                    # Standardising Research Summary
                    std_research_summary_df = standardise_and_transform_research(dataframes["Research Summary"])
                    std_research_summary_df = std_research_summary_df.loc[:, std_research_summary_df.columns.intersection(standard_columns)]
                    std_research_summary_df = std_research_summary_df.reindex(columns=standard_columns)


                    # Standardising ATM Summary
                    std_atm_summary_df = standardise_and_transform_atm(dataframes["ATM Summary"])
                    std_atm_summary_df = std_atm_summary_df.loc[:, std_atm_summary_df.columns.intersection(standard_columns)]
                    std_atm_summary_df = std_atm_summary_df.reindex(columns=standard_columns)


                    # Combine all standardised dataframes
                    combined_df = pd.concat([
                        std_lab_cases_df, std_invitae_summary_df, std_clinical_summary_df, std_research_summary_df, std_atm_summary_df
                    ], axis=0, ignore_index=True)
                    st.session_state["Combined"] = combined_df
                    st.success("✅ Data successfully standardised and combined!")
                    
                    # Display results
                    st.write("### Combined Dataframe")
                    st.write(combined_df)
            except Exception as e:
                    st.error(f"❌ Failed to standardise data: {e}")
        else:
            st.warning("⚠️ Please load data first.")



    # Query Data Tab
    elif variant_tab == "Query Data":
        st.header("🔍 Query Variant Data")

        if "Combined" in st.session_state:
            query(st.session_state["Combined"])
        else:
            st.warning("⚠️ Please standardise and combine data before querying.")



    # Validate Variants Tab
    elif variant_tab == "Validate Variants":
        st.header("✅ Variant Validator")

        uploaded_file = st.file_uploader("Upload CSV file with variants", type=["csv"])
        if uploaded_file:
            df = pd.read_csv(uploaded_file)
            if 'Transcript' not in df.columns or 'Variant (HGVSc)' not in df.columns:
                st.error("❌ Uploaded file must contain 'Transcript' and 'Variant (HGVSc)' columns.")
            else:
                st.success("✅ File uploaded successfully!")
                df['HGVS (HGVSc)'] = df['Transcript'] + ":" + df['Variant (HGVSc)']

                if st.button("Validate Variants"):
                    with st.spinner("Validating variants..."):
                        results_df = test_variant_validator_batch(df['HGVS (HGVSc)'].tolist())

                    st.success("✅ Validation completed!")
                    st.dataframe(results_df)
                    csv_data = results_df.to_csv(index=False)
                    st.download_button("⬇️ Download Validation Results (CSV)", csv_data, "validation_results.csv", "text/csv")

# Footer
    st.markdown("---")
    st.markdown("💡 **Tip:** Upload `config.yaml`, standardise the data, then query results!")

# Tab 2: FastHPOCR
with tab2: 
    # Sidebar: Upload HPO Ontology File
    st.sidebar.subheader("HPO Annotation")
    uploaded_file = st.sidebar.file_uploader("Upload `hp.obo` file", type=["obo"])

    if uploaded_file:
        # Save the uploaded file
        hpo_path = os.path.join(OUTPUT_FOLDER, "hp.obo")
        os.makedirs(OUTPUT_FOLDER, exist_ok=True)

        with open(hpo_path, "wb") as f:
            f.write(uploaded_file.read())

        st.sidebar.success(f"`hp.obo` uploaded successfully!")

        # Indexing Configuration
        st.sidebar.header("⚙️ Indexing Options")

        # Root Concept Selection (Multi-select dropdown)
        root_concept_choices = {
            "HP:0000118": "Phenotypic abnormality",
            "HP:0000707": "Nervous system",
            "HP:0001574": "Growth",
            "HP:0000478": "Eye",
            "HP:0000598": "Ear",
            "HP:0000119": "Metabolism/Homeostasis",
            "HP:0003011": "Muscle",
            "HP:0000769": "Voice",
            "HP:0001507": "Skeletal",
            "HP:0001939": "Blood & Hematology",
            "HP:0001626": "Cardiovascular",
            "HP:0000818": "Breast",
            "HP:0001392": "Digestive system",
            "HP:0002086": "Respiratory",
            "HP:0002715": "Endocrine",
            "HP:0000924": "Genitourinary",
            "HP:0002664": "Immune system",
            "HP:0000767": "Integument (Skin, Hair, Nails)",
            "HP:0025031": "Prenatal & Birth",
            "HP:0001197": "Neoplasm",
            "HP:0032443": "Cellular phenotype",
            "HP:0001871": "Immune system",
            "HP:0000152": "Mouth & Dentition",
            "HP:0011446": "Connective tissue"
        }

        selected_root_concepts = st.sidebar.multiselect(
            "Select Root Concepts",
            options = list(root_concept_choices.keys()),
            default = list(root_concept_choices.keys())  # Pre-select all
        )

        # Boolean Configuration Settings
        allow_acronyms = st.sidebar.checkbox("Allow 3-letter acronyms (e.g., 'VUR')", value=True)
        include_categories = st.sidebar.checkbox("Include top-level category", value=True)
        allow_duplicates = st.sidebar.checkbox("Allow duplicate entries", value=False)

        if st.sidebar.button("🔧 Generate Index"):
            with st.spinner("Indexing HPO terms..."):
                indexConfig = {
                    "rootConcepts": selected_root_concepts,
                    "allow3LetterAcronyms": allow_acronyms,
                    "includeTopLevelCategory": include_categories,
                    "allowDuplicateEntries": allow_duplicates
                }

                indexHPO = IndexHPO(hpo_path, OUTPUT_FOLDER, indexConfig=indexConfig)
                indexHPO.index()

            st.sidebar.success("✅ HPO Index Created!")

    # Main Section: Text Annotation
    st.header("📝 Annotate Medical Text")
    text_input = st.text_area("Enter medical text", "The patient was diagnosed with vesicoureteral reflux and muscle weakness.")

    if st.button("🔍 Annotate"):
        if not os.path.exists(INDEX_FILE):
            st.error("❌ No index found! Please upload `hp.obo` and generate the index first.")
        else:
            with st.spinner("Processing..."):
                annotator = HPOAnnotator(INDEX_FILE)
                annotations = annotator.annotate(text_input)

            if annotations:
                st.success("✅ HPO Concepts Found!")

                # Convert annotations to structured data
                annotation_data = []
                for ann in annotations:
                    annotation_data.append({
                        "HPO ID": ann.hpoUri,
                        "Label": ann.hpoLabel,
                        "Text Span": ann.textSpan
                    })
                    st.write(f"🧬 **{ann.hpoUri}** - {ann.hpoLabel} (Found: '{ann.textSpan}')")

                # Convert to Dataframe
                df = pd.DataFrame(annotation_data)

                # JSON Download Button
                json_filename = "annotations.json"
                with open(json_filename, "w") as f:
                    json.dump(annotation_data, f)

                st.download_button(
                    label="⬇️ Download Annotations (JSON)",
                    data=json.dumps(annotation_data, indent=4),
                    file_name=json_filename,
                    mime="application/json"
                )

                # CSV Download Button
                csv_filename = "annotations.csv"
                df.to_csv(csv_filename, index=False)

                st.download_button(
                    label="⬇️ Download Annotations (CSV)",
                    data=df.to_csv(index=False).encode("utf-8"),
                    file_name=csv_filename,
                    mime="text/csv"
                )

            else:
                st.warning("⚠️ No HPO concepts detected.")

    # Footer
    st.markdown("---")
    st.markdown("💡 **Tip:** Upload `hp.obo`, generate the index, then annotate text!")


# **Tab 3: About**
with tab3:
    st.header("📚 About")
    st.markdown("""
    The **PhenoVariant App** is a modular tool for querying, annotating, and validating
    genomic variants while ensuring privacy, preserving data access.

    The **Variant Database** was developed by the 2025 New Colombo Plan interns, Angelo Lagahit and Sze Wei Shong, from Curtin University.

    The **HPO Annotation** model was developed by Tudor Groza, and integrated by the interns.
    """, unsafe_allow_html=True)

    st.subheader("🧬 Variant Database")
    st.markdown("""
    - Consolidates & harmonises variant data  
    - Supports multi-level querying (Gene + Phenotype + Variant Type + Solved Status)  
    - Layered filtering, export, and external validation (VariantValidator)""")

    st.subheader("🔍 HPO Text Annotator")
    st.markdown("""
    - Extracts HPO concepts from free text with high-speed NLP (FastHPOCR)""")

