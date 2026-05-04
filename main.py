"""Entry point for Streamlit Cloud — delegates to app/main.py."""                                                                                       
import os                                                 
import sys

# app/main.py uses unprefixed imports (`from core.checker import ...`),                                                                                 
# expecting `app/` on sys.path. Mirror what `streamlit run app/main.py`
# would do so the shim behaves identically.                                                                                                             
HERE = os.path.dirname(os.path.abspath(__file__))         
sys.path.insert(0, os.path.join(HERE, "app"))     

import app.main  # noqa: F401
