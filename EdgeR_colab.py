import subprocess

import pandas as pd
import ipywidgets as widgets
from IPython.display import display

counts_path = widgets.Text(description="Counts CSV:")
metadata_path = widgets.Text(description="Metadata CSV:")

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
            "EdgeR_pipeline.R",
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
