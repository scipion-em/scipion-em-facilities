# quality_dashboard.py
import os
import streamlit as st
import pandas as pd
from streamlit_autorefresh import st_autorefresh
import altair as alt
import matplotlib.pyplot as plt
from collections import OrderedDict
import mrcfile
import numpy as np
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import re
from st_aggrid.shared import JsCode


VIEWS = ('Main View', 'Dose', 'Drift', 'Tilt', 'Miffi', 'CTF', 'Micrographs scores')
CREATION_TIME = 'creationTime'
ACQUISITION_COLUMS = ['magnification','pixelSize','voltage','sphericalAberration','dosePerFrame']

ALL_COLUMNS = {
    'import': ['movieId','movieName','creationTime','magnification','pixelSize','voltage','sphericalAberration','dosePerFrame'],
    'dose': ['movieId','movieName','Filter DoseAnalysis','thresholdPercentageDiff','diffDosePerAngstrom2','meanDosePerAngstrom2','stdDosePerAngstrom2'],
    'drift': ['movieId','movieName','Filter MaxShift','thresholdMaxMovieShift','maxMovieShift','thresholdMaxFrameShift','maxFrameShift','accumMotionTotal','accumMotionEarly','accumMotionLate','micName'],
    'tilt': ['movieId','movieName','Filter TiltAnalysis','thresholdMeanCorrelation','tiltMeanCorrelation','thresholdStdCorrelation','tiltStdCorrelation','tiltImage','micName'],
    'miffi': ['movieId','movieName','Filter Miffi','miffiLabel','micName'],
    'ctf':['movieId','movieName','Filter CTFConsensus','thresholdResolution','resolution','fitQuality','thresholdAstigmatismPercentage','astigmatismPercentage','defocusU','defocusV','defocusRatio','defocusAngle','IceRingDensity','thresholdMaxDefocus','thresholdMinDefocus','consensusResolution','thresholdConsensusResolution','micName','psdFile']
}

VISIBLE_COLUMNS = {
    'import': ['movieId','movieName','creationTime'],
    'dose': ['movieId','movieName','creationTime','Filter DoseAnalysis','diffDosePerAngstrom2','meanDosePerAngstrom2','stdDosePerAngstrom2'],
    'drift': ['movieId','micName','creationTime','Filter MaxShift','maxFrameShift','maxMovieShift','accumMotionTotal','accumMotionEarly','accumMotionLate'],
    'tilt': ['movieId','micName','creationTime','Filter TiltAnalysis','tiltMeanCorrelation','tiltStdCorrelation'],
    'miffi': ['movieId','micName','creationTime','Filter Miffi','miffiLabel'],
    'ctf':['movieId','micName','creationTime','Filter CTFConsensus','resolution','fitQuality', 'defocusRatio', 'defocusU','defocusV','astigmatismPercentage','defocusAngle','IceRingDensity','consensusResolution']
}

def load_data(path):
    try:
        df = pd.read_csv(path)
        # Make feature engineering
        df[CREATION_TIME] = pd.to_datetime(df[CREATION_TIME])
        # Rename 'boolSomething' to 'Filter Something'
        df = df.rename(columns={
            col: 'Filter ' + col[8:]
            for col in df.columns if col.startswith('boolPass')
        })
        # Round all float columns to 2 decimals
        df = df.round(2)
        return df
    except Exception as e:
        st.error(f"Error loading CSV: {e}")
        return pd.DataFrame()

def plot_histogram_with_threshold(df, column, threshold_col_min=None, threshold_col_max=None, bins=50):
    if column not in df.columns:
        return
    st.markdown(f"#### {column}")
    hist_values = df[column].dropna()
    fig = None
    try:
        fig, ax = plt.subplots()
        ax.hist(hist_values, bins=bins, color='skyblue', edgecolor='black')
        if threshold_col_min and threshold_col_min in df.columns:
            ax.axvline(df[threshold_col_min].iloc[0], color='red', linestyle='--', label='Min Threshold')
        if threshold_col_max and threshold_col_max in df.columns:
            ax.axvline(df[threshold_col_max].iloc[0], color='green', linestyle='--', label='Max Threshold')
        ax.set_xlabel(column)
        ax.set_ylabel("Count")
        ax.legend()
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error plotting {column}: {e}")

def plot_filters_distribution(df, label, title):
    title_hist = title + ' Distribution'
    # st.subheader(title_hist)
    st.markdown(f"##### {title_hist}")

    # Count accepted and rejected
    filters_counts = df[label].value_counts().rename({True: 'Accepted', False: 'Rejected'})
    total = filters_counts.sum()

    # Create a DataFrame for plotting
    plot_df = filters_counts.reset_index()
    plot_df.columns = ['Status', 'Count']
    plot_df['Percentage'] = (plot_df['Count'] / total * 100).round(1)
    plot_df['Color'] = plot_df['Status'].map({'Accepted': 'green', 'Rejected': 'red'})

    # Create Altair chart
    chart = alt.Chart(plot_df).mark_bar().encode(
        x=alt.X('Status:N', title=""),
        y=alt.Y('Count:Q', title="Number of Images"),
        color=alt.Color('Status:N', scale=alt.Scale(domain=['Accepted', 'Rejected'], range=['green', 'red']),
                        legend=None)
    ).properties(width=300, height=500)

    # Add percentage labels
    text = alt.Chart(plot_df).mark_text(
        align='center',
        baseline='bottom',
        dy=-5
    ).encode(
        x='Status:N',
        y='Count:Q',
        text=alt.Text('Percentage:Q', format='.1f')
    )

    st.altair_chart(chart + text, use_container_width=False)

def generate_js_threshold_vars(cond):
    # Extract and declare JS variables from threshold placeholders
    matches = re.findall(r'threshold\w+', cond)
    lines = []
    for var in matches:
        lines.append(f"let {var} = row['{var}'];")
    return '\n'.join(lines)

def replace_threshold_vars(cond):
    return cond.replace("value", "value")

def configure_aggrid_threshold_styles(gb, thresholds_conditions):
    for col, condition_str in thresholds_conditions.items():

        # --- Miffi filter condition ---
        if condition_str.startswith("Filter Miffi"):
            gb.configure_column(
                col,
                cellStyle=JsCode("""
                    function(params) {
                        if (params.value === true) {
                            return {backgroundColor: 'lightgreen', color: 'black'};
                        } else if (params.value === false) {
                            return {backgroundColor: 'lightcoral', color: 'white'};
                        } else {
                            return {};
                        }
                    }
                """)
            )
            continue
        # --- Numeric threshold logic ---
        js_vars = generate_js_threshold_vars(condition_str)
        js_condition = replace_threshold_vars(condition_str)
        js_code = f"""
            function(params) {{
                const row = params.data;
                let value = row['{col}'];
                if (value == null || isNaN(value)) return;
                {js_vars}
                if ({js_condition}) {{
                    return {{backgroundColor: 'lightgreen', color: 'black'}};
                }} else {{
                    return {{backgroundColor: 'lightcoral', color: 'white'}};
                }}
            }}
        """
        gb.configure_column(col, cellStyle=JsCode(js_code))

def main_view(df):
    st.subheader("Filters statistics")
    bool_cols = [col for col in df.columns if col.startswith("Filter")]
    # Group into rows of 2 columns
    for i in range(0, len(bool_cols), 2):
        cols_row = bool_cols[i:i + 2]
        columns_st = st.columns(len(cols_row))
        for col_st, col in zip(columns_st, cols_row):
            with col_st:
                passed = df[col].sum()
                total = df[col].count()
                percentage = (passed / total) * 100 if total > 0 else 0

                color = "green" if percentage >= 60 else "red"

                st.markdown(f"""
                    <div style="text-align: center; padding: 10px; border: 1px solid #eee; border-radius: 5px;">
                        <strong>{col}</strong><br>
                        <span style="color: {color}; font-size: 18px;">{passed}/{total} ({percentage:.1f}%)</span>
                    </div>
                """, unsafe_allow_html=True)

    st.subheader("General dataframe")
    # Flatten and remove duplicates while preserving order
    all_visible_columns = list(OrderedDict.fromkeys(
        col for cols in VISIBLE_COLUMNS.values() for col in cols))

    df_filter_visible = df[[col for col in all_visible_columns if col in df.columns]].copy()
    st.dataframe(df_filter_visible, use_container_width=True)

def common_view_filters(df, label, title, columns_tag, thresholds_conditions):
    # --------------------  DEFINE WINDOW STRUCTURE  --------------------------------------
    # Top Panel
    with st.container():
        # Split into two horizontal panels (left and right)
        container_top_left, container_top_right = st.columns(2)
        with container_top_left:
            panel_distribution, panel_thresholds_button = st.columns(2)

    # Bottom panel
    with st.container():
        tab_table, tab_plots = st.tabs(["Dataframe", "Plots"])

    if label not in df.columns:
        st.write(f'No views for the {title} yet')
        return

    with panel_distribution:
        plot_filters_distribution(df, label, title)

    with panel_thresholds_button:
        with st.container():
            threshold_values = {}
            miffi_labels = []
            for col, condition in thresholds_conditions.items():
                if condition.startswith("Filter Miffi"):
                    filter_col = 'miffiLabel'
                    if filter_col in df.columns:
                        miffi_labels = df[filter_col].dropna().unique().tolist()
                else:
                    threshold_vars = [var for var in df.columns if var.startswith("threshold") and var in condition]
                    for var in threshold_vars:
                        if var not in threshold_values:
                            if var in df.columns and not df[var].dropna().empty:
                                threshold_values[var] = df[var].dropna().iloc[0]
                            else:
                                threshold_values[var] = "N/A"

            # Display Miffi labels if any
            if miffi_labels:
                st.markdown("##### Miffi Labels Present in Dataset:")
                st.write(", ".join(map(str, miffi_labels)))
            else:
                st.markdown("##### Thresholds used for filtering:")
                for idx, (col, val) in enumerate(threshold_values.items()):
                    with st.container():
                        st.metric(label=col, value=f"{val}")

        with st.container():
            # Interactive control to filter
            filter_option = st.radio("Show data for:", ["All", "Accepted", "Rejected"])

            if filter_option == "Accepted":
                df = df[df[label] == True]
            elif filter_option == "Rejected":
                df =  df[df[label] == False]
            else:
                df = df

    df_filter_all = df[[col for col in ALL_COLUMNS[columns_tag] if col in df.columns]].copy()

    visible_cols = [col for col in VISIBLE_COLUMNS[columns_tag] if col in df.columns]
    all_cols = df_filter_all.columns.tolist()
    hidden_cols = [col for col in all_cols if col not in visible_cols]

    with tab_table:
        # Build interactive grid
        gb = GridOptionsBuilder.from_dataframe(df_filter_all)
        gb.configure_selection('single')  # one row selectable at a time
        # Apply conditional formatting styles
        configure_aggrid_threshold_styles(gb, thresholds_conditions)

        # Hide all columns not in VISIBLE_COLUMNS
        for col in hidden_cols:
            gb.configure_column(col, hide=True)
        # Create the grid
        grid_options = gb.build()
        # Display interactive grid
        grid_response = AgGrid(
            df_filter_all,
            gridOptions=grid_options,
            update_mode=GridUpdateMode.SELECTION_CHANGED,
            allow_unsafe_jscode=True,
            enable_enterprise_modules=False,
            height=400,
            theme="streamlit",  # "material", "alpine", etc.
        )
        # Get selected row (if any)
        mic_path= ''
        selected_rows = pd.DataFrame(grid_response['selected_rows'])

        if not selected_rows.empty:
            mic_path = selected_rows.iloc[0]["micName"]
            st.success(f"Selected: {mic_path}")

    with container_top_right:
        st.markdown("##### Micrograph Viewer")
        # Let user select a micrograph to view
        if os.path.exists(mic_path):
            with mrcfile.open(mic_path, permissive=True) as mrc:
                data = mrc.data
            # Normalize if necessary
            image = (data - np.min(data)) / (np.max(data) - np.min(data))  # scale 0-1
            st.image(image, caption=os.path.basename(mic_path), clamp=True)
        else:
            print('Not an image selected')

    return df, tab_plots

def dose_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'diffDosePerAngstrom2': "value < thresholdPercentageDiff",
    }
    df_filter, tab_plots = common_view_filters(df, label, title, columns_tag, thresholds_conditions)

def drift_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'maxMovieShift': "value < thresholdMaxMovieShift",
        'maxFrameShift': "value < thresholdMaxFrameShift",
    }
    df_filter, tab_plots = common_view_filters(df, label, title, columns_tag, thresholds_conditions)

def tilt_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'tiltMeanCorrelation': "thresholdMeanCorrelation < value",
        'tiltStdCorrelation': "value < thresholdStdCorrelation",
    }
    df_filter, tab_plots = common_view_filters(df, label, title, columns_tag, thresholds_conditions)

def miffi_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'miffiLabel': "Filter Miffi"
    }
    df_filter, tab_plots = common_view_filters(df, label, title, columns_tag, thresholds_conditions)


def ctf_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'defocusU': "thresholdMinDefocus < value < thresholdMaxDefocus",
        'defocusV': "thresholdMinDefocus < value < thresholdMaxDefocus",
        'astigmatismPercentage': "value < thresholdAstigmatismPercentage",
        'resolution': "value < thresholdResolution",
        'consensusResolution': "value < thresholdConsensusResolution",
    }

    df, tab_plots = common_view_filters(df, label, title, columns_tag, thresholds_conditions)

    with tab_plots:
        plot_histogram_with_threshold(df,"defocusU", threshold_col_min="thresholdMinDefocus",
                                      threshold_col_max="thresholdMaxDefocus")
        plot_histogram_with_threshold(df,"defocusV", threshold_col_min="thresholdMinDefocus",
                                      threshold_col_max="thresholdMaxDefocus")
        plot_histogram_with_threshold(df,"astigmatismPercentage", threshold_col_max="thresholdAstigmatismPercentage")
        plot_histogram_with_threshold(df,"resolution", threshold_col_max="thresholdResolution")
        plot_histogram_with_threshold(df,"consensusResolution", threshold_col_max="thresholdConsensusResolution")

        # Scatter plot DefocusU vs DefocusV
        if 'defocusU' in df.columns and 'defocusV' in df.columns:
            st.markdown("### DefocusU vs DefocusV (Astigmatism Check)")
            try:
                import altair as alt
                scatter_df = df[['defocusU', 'defocusV']].dropna()
                chart = alt.Chart(scatter_df).mark_circle(size=60, opacity=0.5).encode(
                    x='defocusU:Q',
                    y='defocusV:Q',
                    tooltip=['defocusU', 'defocusV']
                ).properties(
                    width=600,
                    height=400
                )
                st.altair_chart(chart, use_container_width=True)
            except Exception as e:
                st.error(f"Error plotting defocus scatter: {e}")

def micrograph_view(df):
    # Sidebar
    st.sidebar.title("Settings")
    selected_label = st.sidebar.selectbox("Choose filter label", df.columns)  # Very interesting for selecting one column as label

    # Horizontal panels
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Filter Distribution")
        plot_filters_distribution(df, selected_label, "Filtering Status")

    with col2:
        threshold_values = {
            'threshold1':1,
            'threshold2':2
        }
        st.subheader("Thresholds Summary")
        with st.container():
            for threshold, value in threshold_values.items():
                st.metric(label=threshold, value=value)

    # Tabs for detailed views
    tab1, tab2 = st.tabs(["Filtered View", "All Data"])

    with tab1:
        st.subheader("Filtered Data")
        # common_view_filters(df[df[selected_label] == True], selected_label, "Accepted", "your_tag",
        #                     thresholds_conditions)

    with tab2:
        st.subheader("All Data")
        # common_view_filters(df, selected_label, "All", "your_tag", thresholds_conditions)

    # Expanders for optional content
    with st.expander("Show Debug Info"):
        st.write(df.describe())

def main():
    st.set_page_config(page_title="Quality Monitor", layout="wide")
    st.title("Live Quality Metrics Monitor")

    # Auto-refresh every 30 seconds (30,000 milliseconds)
    st_autorefresh(interval=30_000, limit=None, key="quality_monitor_refresh")

    monitor_path = os.environ.get("STREAMLIT_MONITOR_PATH") # .csv metadata

    if not monitor_path or not os.path.exists(monitor_path):
        st.error("The CSV metadata file does not exists.")
        return

    df = load_data(monitor_path)

    if df.empty:
        st.warning("The CSV is empty or it could not be loaded.")
        return

    # Select the last 10 minutes
    df_diff = df.copy()
    metric_date_10min = df_diff['creationTime'].max() - pd.DateOffset(minute=10)

    # --------------  BUILD DASHBOARD ------------------------------
    add_sidebar = st.sidebar.selectbox('Objects Views', VIEWS)
    # ('Main View', 'Dose', 'Drift', 'Tilt', 'Miffi', 'CTF', 'Micrographs scores')
    if add_sidebar == 'Main View':
        main_view(df)

    if add_sidebar == 'Dose':
        dose_view(df,'Filter DoseAnalysis', 'Dose Analysis Filter', 'dose')

    if add_sidebar == 'Drift':
        drift_view(df,'Filter MaxShift', 'Drift Analysis Filter', 'drift')

    if add_sidebar == 'Tilt':
        tilt_view(df, 'Filter TiltAnalysis', 'Tilt Analysis Filter', 'tilt')

    if add_sidebar == 'Miffi':
        miffi_view(df, 'Filter Miffi', 'Miffi Filter', 'miffi')

    if add_sidebar == 'CTF':
        ctf_view(df, 'Filter CTFConsensus', 'CTF Filter', 'ctf')

    if add_sidebar == 'Micrographs scores':
        micrograph_view(df)

if __name__ == "__main__":
    main()
