import subprocess
import os

import pandas as pd
import ipywidgets as widgets
from IPython.display import display
from IPython import get_ipython

user_ns = get_ipython().user_ns if get_ipython() else {}
default_counts = user_ns.get("EXP_FILE") or os.environ.get("EXP_FILE", "")
default_metadata = user_ns.get("META_FILE") or os.environ.get("META_FILE", "")

counts_path = widgets.Text(description="Counts CSV:", value=default_counts)
metadata_path = widgets.Text(description="Metadata CSV:", value=default_metadata)

sample_dropdown = widgets.Dropdown(description="Sample column:")
group_dropdown = widgets.Dropdown(description="Group column:")


def update_columns(change):
    path = metadata_path.value
    if path:
        try:
            df = pd.read_csv(path, nrows=0)
            cols = df.columns.tolist()
            sample_dropdown.options = cols
            group_dropdown.options = cols
        except Exception as e:
            sample_dropdown.options = []
            group_dropdown.options = []


metadata_path.observe(update_columns, names="value")

update_columns(None)

run_button = widgets.Button(description="Run EdgeR pipeline")
output = widgets.Output()


def run_pipeline(_):
    with output:
        output.clear_output()
        if not (
            counts_path.value
            and metadata_path.value
            and sample_dropdown.value
            and group_dropdown.value
        ):
            print("Please provide all inputs before running the pipeline.")
            return

        cmd = [
            "Rscript",
            "/content/EdgeR_pipeline.R",
            counts_path.value,
            metadata_path.value,
            sample_dropdown.value,
            group_dropdown.value,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)


run_button.on_click(run_pipeline)

display(counts_path, metadata_path, sample_dropdown, group_dropdown, run_button, output)
