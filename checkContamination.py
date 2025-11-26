import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from matplotlib_venn import venn2, venn3


class ContaminationChecker:
    """Class to perform contamination checks and filtering on fungi data."""

    def __init__(self, curated_df, score_weights):
        self.curated_df = curated_df
        self.default_score_weights = score_weights

    def flatten_set_of_lists(self, set_of_lists):
        flattened_list = [item for sublist in set_of_lists for item in sublist]
        return set(flattened_list)

    def get_unique_properties(self, filtered_fungi):
        # Get unique properties from the 'Contributing Properties' column
        all_properties = self.flatten_set_of_lists(
            filtered_fungi["Contributing Properties"].dropna()
        )
        if not all_properties:
            print("Debug: all_properties set is empty")
            print("Sample of 'Contributing Properties' column:")
            print(filtered_fungi["Contributing Properties"].head())
        return list(all_properties)

    def normalize_species_name(self, name):
        """Normalize species name by removing sp., spp. and extracting genus if needed."""
        if pd.isna(name):
            return None, False
        name = name.strip()
        lower = name.lower()
        if lower.endswith(" sp.") or lower.endswith(" spp.") or lower.endswith(" sp") or lower.endswith(" spp"):
            genus = name.split()[0]
            return genus, True  # group
        else:
            return name, False  # exact

    def generate_venn_diagram(self, filtered_fungi):
        """Generate a Venn diagram of contributing properties."""
        properties = self.get_unique_properties(filtered_fungi)
        num_properties = len(properties)

        if num_properties == 0:
            st.write("No contributing properties found.")
            return None
        elif num_properties == 1:
            count = sum(
                properties[0] in item
                for item in filtered_fungi["Contributing Properties"]
                if isinstance(item, list)
            )
            st.write(
                f"Only one property found: {properties[0]}. {count} Fungi belong to this property."
            )
            return None
        elif num_properties == 2:
            subsets = [
                set(
                    filtered_fungi[
                        filtered_fungi["Contributing Properties"].apply(
                            lambda x: prop in x
                        )
                    ].index
                )
                for prop in properties
            ]
            plt.figure(figsize=(4, 4))
            venn2(subsets=subsets, set_labels=properties)
        else:
            col1, col2, col3 = st.columns(3)
            prop1 = col1.selectbox("Property 1", properties, index=0)
            prop2 = col2.selectbox("Property 2", properties, index=1)
            prop3 = col3.selectbox(
                "Property 3",
                ["None"] + properties,
                index=3 if num_properties > 2 else 0,
            )

            if len(set([prop1, prop2, prop3])) < 3 or prop3 == "None":
                if prop1 == prop2 or (
                    prop3 != "None" and (prop1 == prop3 or prop2 == prop3)
                ):
                    st.warning("Please select distinct properties.")
                    return None
                subsets = [
                    set(
                        filtered_fungi[
                            filtered_fungi["Contributing Properties"].apply(
                                lambda x: prop in x
                            )
                        ].index
                    )
                    for prop in [prop1, prop2]
                ]
                plt.figure(figsize=(4, 4))
                venn2(subsets=subsets, set_labels=[prop1, prop2])
            else:
                subsets = [
                    set(
                        filtered_fungi[
                            filtered_fungi["Contributing Properties"].apply(
                                lambda x: prop in x
                            )
                        ].index
                    )
                    for prop in [prop1, prop2, prop3]
                ]
                plt.figure(figsize=(4, 4))
                venn3(subsets=subsets, set_labels=[prop1, prop2, prop3])

        plt.title("Contributing Properties of Filtered Fungi")
        return plt.gcf()

    def filter_fungi(
        self,
        input_df,
        score_weights,
        score_threshold,
        reads_threshold,
    ):
        """Filter fungi based on thresholds and weights, with group handling."""
        species_column = input_df.columns[0]
        self.species_column_name = species_column
        curated_species_list = self.curated_df["Species"].tolist()

        # --- 🧠 Group + exact matching logic ---
        group_stats = []
        matched_indices = []
        group_mapping = {}  # map input row -> matched curated species list

        for idx, input_name in input_df[species_column].items():
            norm_name, is_group = self.normalize_species_name(input_name)
            if norm_name is None:
                continue

            if is_group:
                matches = [
                    s for s in curated_species_list
                    if s.lower().startswith(norm_name.lower() + " ")
                ]
                if matches:
                    matched_indices.append(idx)
                    group_mapping[idx] = matches
                    group_stats.append({
                        "Group": input_name,
                        "Genus": norm_name,
                        "Number of matches": len(matches),
                        "Matched curated species": matches
                    })
            else:
                if input_name in curated_species_list:
                    matched_indices.append(idx)
                    group_mapping[idx] = [input_name]

        matching_rows_df = input_df.loc[matched_indices]
        # st.dataframe(matching_rows_df)
        unmatched_indices = list(set(input_df.index) - set(matched_indices))
        non_matching_rows_df = input_df.loc[unmatched_indices]
        self.non_matching_rows_df = non_matching_rows_df
        self.non_matching_rows = non_matching_rows_df.shape[0]

        # --- 📊 Location filtering ---
        location_columns = input_df.columns[1:]
        matching_rows_df["Num loc"] = matching_rows_df[location_columns].apply(
            lambda x: x[x >= reads_threshold].count(), axis=1
        )
        matching_rows_df["Locations"] = matching_rows_df[location_columns].apply(
            lambda x: {
                loc: count for loc, count in x.items() if count >= reads_threshold
            },
            axis=1,
        )
        matching_rows_df = matching_rows_df[matching_rows_df["Num loc"] > 0]
        if matching_rows_df.empty:
            print("Warning: No Fungi species meet the location count threshold.")
            return 0, 0, 0, 0, pd.DataFrame()
        
        # st.dataframe(matching_rows_df)

        # --- 🧮 Score calculation with average for groups ---
        properties = [p for p in score_weights.keys() if p in self.curated_df.columns]
        weights = [score_weights[p] for p in properties]

        def calc_avg_score(idx):
            matched_species_list = group_mapping[idx]
            scores = []
            contributing_props_all = set()
            for sp in matched_species_list:
                row = self.curated_df[self.curated_df["Species"] == sp]
                if row.empty:
                    continue
                vals = row[properties].values.flatten()
                min_len = min(len(vals), len(weights))
                weighted_vals = vals[:min_len] * weights[:min_len]
                score = sum(weighted_vals)
                scores.append(score)
                props = [prop for prop, val in zip(properties, weighted_vals) if val > 0]
                contributing_props_all.update(props)
            if not scores:
                return 0, []
            avg_score = np.mean(scores)
            return avg_score, list(contributing_props_all)

        scores = []
        props_list = []
        for idx in matching_rows_df.index:
            avg_score, props_used = calc_avg_score(idx)
            scores.append(avg_score)
            props_list.append(props_used)

        matching_rows_df["Weight Score"] = scores
        matching_rows_df["Contributing Properties"] = props_list

        # st.dataframe(matching_rows_df)
        total_matches = len(matching_rows_df)

        # Mark average for groups in output
        def label_group(row):
            norm_name, is_group = self.normalize_species_name(row[species_column])
            if is_group:
                return f"{row['Weight Score']:.2f} (avg)"
            else:
                return f"{row['Weight Score']:.2f}"

        matching_rows_df = matching_rows_df[matching_rows_df["Weight Score"] >= score_threshold]
        if matching_rows_df.empty:
            return 0, 0, 0, 0, pd.DataFrame()

        matching_rows_df["Score"] = matching_rows_df.apply(label_group, axis=1)

        # st.dataframe(matching_rows_df)

        filtered_fungi = matching_rows_df[
            [
                self.species_column_name,
                "Score",
                "Contributing Properties",
                "Num loc",
                "Locations",
            ]
        ]

        # --- Reverse table ---
        unique_props = self.get_unique_properties(filtered_fungi)
        property_species_data = []
        for prop in unique_props:
            matching_species = filtered_fungi[
                filtered_fungi["Contributing Properties"].apply(lambda x: prop in x)
            ][self.species_column_name].tolist()
            property_species_data.append(
                {"Property": prop, "Number": len(matching_species), "Matching Species": matching_species}
            )
        reverse_table = pd.DataFrame(property_species_data)

        # --- 📊 Group statistics table ---
        if group_stats:
            group_stats_df = pd.DataFrame(group_stats)
            # st.subheader("Group Match Statistics")
            group_stats_df = group_stats_df.drop_duplicates(subset=["Group"])
            # st.dataframe(group_stats_df)
        else:
            group_stats_df = pd.DataFrame()

        matching_rows = matching_rows_df.shape[0]
        thresh_rows = matching_rows

        return total_matches, filtered_fungi, thresh_rows, reverse_table, group_stats_df
