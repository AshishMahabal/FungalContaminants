import streamlit as st
import pandas as pd
import json
import io
from checkContamination import ContaminationChecker
from display_utils import display_markdown
from blast_pipeline import run_blast_pipeline

# ---------- CONFIG ----------
st.set_page_config(
    page_title="Fungal Contamination Checker",
    page_icon="🧫",
    layout="wide"
)

st.title("🧫 Fungal Contamination Checker Dashboard")

# ---------- HELPERS ----------

def modify_input(input_df, curated_df):
    col1 = input_df.columns[0]
    if col1.lower() != "#Datasets":
        input_df.rename(columns={col1: "#Datasets"}, inplace=True)

    curated_orgs = curated_df["Species"].to_list()
    input_orgs = input_df["#Datasets"].to_list()

    common_orgs = []
    for org in input_orgs:
        try:
            org = float(org)
        except:
            pass
        if isinstance(org, float):
            continue
        org_split = org.split()
        if len(org_split) >= 2:
            org_name = " ".join(org_split[:2])
            for cur_org in curated_orgs:
                if org_name.lower() in cur_org.lower():
                    common_orgs.append((org, cur_org))
                    break

    replaced_map = {}
    for old, new in common_orgs:
        input_df["#Datasets"] = input_df["#Datasets"].replace(old, new)
        replaced_map[new] = old

    for org in input_orgs:
        try:
            org = float(org)
        except:
            continue
        if isinstance(org, float):
            input_df = input_df[input_df["#Datasets"] != str(org)]

    st.session_state["replaced_map"] = replaced_map

    if len(input_df.columns) == 1:
        # add a sample column named 'sample_loc1' with all values as 1
        input_df["sample_loc1"] = 100
    return input_df

# ---------- LOADERS ----------
@st.cache_data
def load_curated_file(choice):
    file_path = (
        "data/curated_fungi.csv"
        if choice == "ID thresh = 35"
        else "data/curated_fungi_75.csv"
    )
    return pd.read_csv(file_path)

@st.cache_data
def load_score_weights():
    with open("data/score_weights.txt", "r") as f:
        return json.load(f)

# ---------- SESSION ----------
if "score_weights" not in st.session_state:
    st.session_state["score_weights"] = load_score_weights()
if "default_weights" not in st.session_state:
    st.session_state["default_weights"] = load_score_weights()
if "curated_choice" not in st.session_state:
    st.session_state["curated_choice"] = "ID thresh = 35"

# ---------- TABS ----------
tab_input, tab_summary, tab_table, tab_credits = st.tabs([
    "🧾 Input Data", "📊 Summary Results", "📄 Detailed Table", "Credits"
])

# ---------- TAB 1: Input ----------
with tab_input:
    st.header("📋 Input & Configuration")
    subtab1, subtab2, subtab3 = st.tabs(["📂 Input File", "⚖️ Weights", "🎚️ Thresholds"])

    with subtab1:
        st.subheader("Upload Input File")
        default_file = st.checkbox("Use sample input (sample-infile.csv)", True)
        uploaded_file = None
        if not default_file:
            uploaded_file = st.file_uploader("Upload CSV file", type="csv")

        curated_choice = st.radio(
            "Curated Species List",
            ["ID thresh = 35", "ID thresh = 75"],
            index=0,
            horizontal=True
        )
        st.session_state["curated_choice"] = curated_choice
        curated_df = load_curated_file(curated_choice)

    with subtab2:
        st.subheader("⚖️ Contamination Weights")
        new_weights = {}
        for key in st.session_state["score_weights"]:
            new_weights[key] = st.slider(
                f"Weight for {key}", 0, 2, st.session_state["score_weights"][key]
            )
        st.session_state["score_weights"] = new_weights

        restore = st.button("Restore Default Weights")
        if restore:
            st.session_state["score_weights"] = st.session_state["default_weights"].copy()

        custom_weights = st.file_uploader("Upload custom weights JSON", type="json")
        if custom_weights:
            st.session_state["score_weights"] = json.load(custom_weights)

    with subtab3:
        st.subheader("🎚️ Threshold Settings")
        score_threshold = st.slider("Score Threshold", 1, 7, 3)
        reads_threshold = st.slider("Reads Threshold", 1, 10000, 10)

    if st.button("🚀 Run Contamination Analysis"):
        if default_file:
            input_df = pd.read_csv("data/sample-infile.csv")
        elif uploaded_file is not None:
            input_df = pd.read_csv(uploaded_file)
        else:
            st.error("Please upload a CSV file.")
            input_df = None

        if input_df is not None:
            input_df = modify_input(input_df, curated_df)
            checker = ContaminationChecker(curated_df, st.session_state["score_weights"])
            results = checker.filter_fungi(
                input_df,
                st.session_state["score_weights"],
                score_threshold,
                reads_threshold
            )

            st.session_state["contamination_results"] = results
            st.session_state["input_df"] = input_df
            st.session_state["curated_df"] = curated_df
            st.session_state["checker"] = checker
            st.success("Analysis completed! Check results in the next tabs.")

# ---------- TAB 2: Summary ----------
with tab_summary:
    st.header("📊 Summary Results")
    if "contamination_results" in st.session_state:
        if isinstance(st.session_state["contamination_results"][1], int):
            st.warning("No contamination detected based on the provided thresholds.")
        else:
            _, _, thresh_rows, _, group_stats_df = st.session_state["contamination_results"]
            input_df = st.session_state["input_df"]
            col1, col2, col3 = st.columns(3)
            col1.metric("Input Records", len(input_df))
            col2.metric("Above Threshold", thresh_rows)
            col3.metric("Curated Species", len(st.session_state["curated_df"]))

            st.markdown("### 🔬 Optional Visuals")
            show_venn = st.checkbox("Show Venn Diagram of Properties", False)
            if show_venn:
                try:
                    fig = st.session_state["checker"].generate_venn_diagram(
                        st.session_state["contamination_results"][1]
                    )
                    st.pyplot(fig)
                except Exception as e:
                    st.error(f"Could not generate Venn: {e}")
            
            if not group_stats_df.empty:
                st.markdown("### 📊 Group Match Statistics")
                group_stats_df = group_stats_df.drop_duplicates(subset=["Group"])
                st.dataframe(group_stats_df)
            else:
                group_stats_df = pd.DataFrame()

    else:
        st.info("Please run the analysis first in the 'Input Data' tab.")

# ---------- TAB 3: Table ----------
with tab_table:
    st.header("📄 Filtered Contamination Table")
    if "contamination_results" in st.session_state:
        if isinstance(st.session_state["contamination_results"][1], int):
            st.warning("No contamination detected based on the provided thresholds.")
        else:
            filtered_df = st.session_state["contamination_results"][1]
            print(filtered_df)
            st.dataframe(filtered_df, use_container_width=True)
            st.download_button(
                "Download CSV",
                data=filtered_df.to_csv(index=False),
                file_name="filtered_fungi.csv",
                mime="text/csv"
            )
    else:
        st.info("Run the analysis first.")

# ---------- TAB 4: Credits ----------
with tab_credits:
    # st.header("Credits")
    try:
        # with open("CREDITS.md", "r", encoding="utf-8") as f:
        #     credits_md = f.read()
        # Use your custom markdown display helper
        display_markdown("CREDITS.md")
        # or, if you prefer plain Streamlit:
        # st.markdown(credits_md)
    except FileNotFoundError:
        st.error("credits.md file not found in the app directory.")

# # ---------- TAB 5: BLAST ----------
# with tab_blast:
#     st.header("🧬 BLAST Pipeline")
#     if "contamination_results" in st.session_state:
#         if isinstance(st.session_state["contamination_results"][1], int):
#             st.warning("No contamination detected; BLAST not applicable.")
#
#         else:
#             checker = st.session_state["checker"]
#             df_unmatched = checker.non_matching_rows_df
#             st.write("Top unmatched species:")
#             st.dataframe(df_unmatched.head(10))
#             if st.button("Run BLAST on Unmatched"):
#                 with st.spinner("Running BLAST..."):
#                     results = run_blast_pipeline(df_unmatched.iloc[:, 0].to_list())
#                 st.success("BLAST complete.")
#     else:
#         st.info("Please run analysis first.")
