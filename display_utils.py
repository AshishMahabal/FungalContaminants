import streamlit as st


def display_markdown(file_path):
    """Read and display a Markdown file."""
    with open(file_path, "r") as file:
        content = file.read()
        st.markdown(content)
