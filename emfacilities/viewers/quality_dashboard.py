# quality_dashboard.py
import os
import streamlit as st
import pandas as pd
from streamlit_autorefresh import st_autorefresh
import altair as alt
import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import OrderedDict
import mrcfile
import numpy as np
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import re
from st_aggrid.shared import JsCode
from streamlit_plotly_events import plotly_events
import plotly.graph_objects as go
from typing import Callable
import seaborn as sns
from itertools import chain


VIEWS = ('Main View', 'Dose', 'Drift', 'Tilt', 'Miffi', 'CTF', 'Scores View')
CREATION_TIME = 'creationTime'
ACQUISITION_COLUMS = ['magnification','pixelSize','voltage','sphericalAberration','dosePerFrame']

ALL_COLUMNS = {
    'import': ['movieId','movieName','creationTime','magnification','pixelSize','voltage','sphericalAberration','dosePerFrame'],
    'dose': ['movieId','micName','creationTime','Filter DoseAnalysis','thresholdPercentageDiff','diffDosePerAngstrom2','meanDosePerAngstrom2','stdDosePerAngstrom2'],
    'drift': ['movieId','micName','creationTime','Filter MaxShift','thresholdMaxMovieShift','maxMovieShift','thresholdMaxFrameShift','maxFrameShift','accumMotionTotal','accumMotionEarly','accumMotionLate'],
    'tilt': ['movieId','micName','creationTime','Filter TiltAnalysis','thresholdMeanCorrelation','tiltMeanCorrelation','thresholdStdCorrelation','tiltStdCorrelation','tiltImage'],
    'miffi': ['movieId','micName','creationTime','Filter Miffi','miffiLabel'],
    'ctf':['movieId','micName','creationTime','Filter CTFConsensus','thresholdResolution','resolution','fitQuality','thresholdAstigmatismPercentage','astigmatismPercentage','defocusU','defocusV','defocusRatio','defocusAngle','IceRingDensity','thresholdMaxDefocus','thresholdMinDefocus','consensusResolution','thresholdConsensusResolution','psdFile']
}

VISIBLE_COLUMNS = {
    'import': ['movieId','movieName','creationTime'],
    'dose': ['movieId','micName','creationTime','Filter DoseAnalysis','diffDosePerAngstrom2','meanDosePerAngstrom2','stdDosePerAngstrom2'],
    'drift': ['movieId','micName','creationTime','Filter MaxShift','maxFrameShift','maxMovieShift','accumMotionTotal','accumMotionEarly','accumMotionLate'],
    'tilt': ['movieId','micName','creationTime','Filter TiltAnalysis','tiltMeanCorrelation','tiltStdCorrelation'],
    'miffi': ['movieId','micName','creationTime','Filter Miffi','miffiLabel'],
    'ctf':['movieId','micName','creationTime','Filter CTFConsensus','resolution','fitQuality', 'defocusRatio', 'defocusU','defocusV','astigmatismPercentage','defocusAngle','IceRingDensity','consensusResolution']
}


# ---------------------------------------- STYLERS -------------------------------------------------------
def generate_js_threshold_vars(condition_str, col):
    matches = re.findall(r'(-)?(threshold\w+)', condition_str)
    js_lines = []
    for sign, var in matches:
        js_var = f"{'neg_' if sign else ''}{var}_{col}"
        accessor = f"row['{var}']"
        if sign:
            js_lines.append(f"let {js_var} = -1 * {accessor};")
        else:
            js_lines.append(f"let {js_var} = {accessor};")
    return "\n".join(js_lines)

def replace_threshold_vars(condition_str, col):
    if '< value <' in condition_str:
        parts = condition_str.split('< value <')
        left = parts[0].strip()
        right = parts[1].strip()
        left_js = re.sub(r'(-)?(threshold\w+)', lambda m: ('neg_' if m.group(1) else '') + m.group(2) + f'_{col}', left)
        right_js = re.sub(r'(-)?(threshold\w+)', lambda m: ('neg_' if m.group(1) else '') + m.group(2) + f'_{col}', right)
        return f"({left_js} < value && value < {right_js})"
    else:
        return re.sub(r'(-)?(threshold\w+)', lambda m: ('neg_' if m.group(1) else '') + m.group(2) + f'_{col}', condition_str)

def configure_aggrid_threshold_styles(gb, thresholds_conditions):
    for col, condition_str in thresholds_conditions.items():
        # --- Miffi filter condition ---
        if condition_str.startswith("Filter Miffi"):
            filter_col = condition_str.strip()
            gb.configure_column(
                col,
                cellStyle=JsCode(f"""
                           function(params) {{
                               let filterVal = params.data["{filter_col}"];
                               if (filterVal === true) {{
                                   return {{backgroundColor: 'lightgreen', color: 'black'}};
                               }} else if (filterVal === false) {{
                                   return {{backgroundColor: 'lightcoral', color: 'white'}};
                               }} else {{
                                   return {{}};
                               }}
                           }}
                       """)
            )
            continue
        # --- Numeric threshold logic ---
        js_vars = generate_js_threshold_vars(condition_str, col)
        js_condition = replace_threshold_vars(condition_str, col)
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

# ---------------------------------------- COMMON FUNCTIONS ----------------------------------------------------
def calculate_filter_stats(df, filter_col):
    passed = df[filter_col].sum()
    total = df[filter_col].count()
    percentage = (passed / total) * 100 if total > 0 else 0
    return passed, total, percentage

def apply_filter_selection(df, label, option):
    if option == "Accepted":
        return df[df[label] == True]
    elif option == "Rejected":
        return df[df[label] == False]
    return df

# ---------------------------------------- MAIN VIEW COMPONENTS ----------------------------------------------------
def render_filter_card(filter_name, passed, total, percentage):
    color = "green" if percentage >= 60 else "red"
    st.markdown(f"""
        <div style="text-align: center; padding: 10px; border: 1px solid #eee; border-radius: 5px;">
            <strong>{filter_name}</strong><br>
            <span style="color: {color}; font-size: 18px;">{passed}/{total} ({percentage:.1f}%)</span>
        </div>
    """, unsafe_allow_html=True)

def render_filter_statistics(df):
    st.subheader("Filters statistics")
    bool_cols = [col for col in df.columns if col.startswith("Filter")]

    for i in range(0, len(bool_cols), 2):
        cols_row = bool_cols[i:i + 2]
        columns_st = st.columns(len(cols_row))
        for col_st, col in zip(columns_st, cols_row):
            with col_st:
                passed, total, percentage = calculate_filter_stats(df, col)
                render_filter_card(col, passed, total, percentage)

def render_one_filter_statistics(df, columns_tag):
    st.markdown("##### Filter statistics")
    df_filtered = df[[col for col in ALL_COLUMNS[columns_tag] if col in df.columns]].copy()
    bool_cols = [col for col in df_filtered.columns if col.startswith("Filter")]
    for col in bool_cols:
        passed, total, percentage = calculate_filter_stats(df_filtered, col)
        render_filter_card(col, passed, total, percentage)
    st.markdown("##### \n")

def render_general_dataframe(df):
    st.subheader("General table")
    all_visible_columns = list(OrderedDict.fromkeys(
        col for cols in VISIBLE_COLUMNS.values() for col in cols))
    df_visible = df[[col for col in all_visible_columns if col in df.columns]].copy()
    st.dataframe(df_visible, use_container_width=True)

# ------------------------ COMMON VIEW STRUCTURE FOR THE DIFFERENT FILTERS ----------------------------------------
def filter_slider_manager(df, visible_vars):
    st.sidebar.markdown("### Interactive plots")
    selected_var = st.sidebar.selectbox("Select variables for interactive plots:", visible_vars)
    df_filtered = df.copy()

    if pd.api.types.is_numeric_dtype(df[selected_var]):
        min_val = float(df[selected_var].min())
        max_val = float(df[selected_var].max())
        step = (max_val - min_val) / 100

        # session_state keys
        min_key = f"{selected_var}_min"
        max_key = f"{selected_var}_max"
        slider_key = f"{selected_var}_slider"

        if st.sidebar.button("Reset to full range"):
            st.session_state[min_key] = min_val
            st.session_state[max_key] = max_val
            st.session_state[slider_key] = (min_val, max_val)

        # Initialized variables
        if min_key not in st.session_state:
            st.session_state[min_key] = min_val
        if max_key not in st.session_state:
            st.session_state[max_key] = max_val

        # Manual inputs
        manual_min = st.sidebar.number_input(
            "Min value",
            value=st.session_state[min_key],
            step=step,
            key=f"{min_key}_input"
        )
        manual_max = st.sidebar.number_input(
            "Max value",
            value=st.session_state[max_key],
            step=step,
            key=f"{max_key}_input"
        )

        # Slider
        selected_range = st.sidebar.slider(
            f"Adjust range of {selected_var}:",
            min_value=min_val,
            max_value=max_val,
            value=(manual_min, manual_max),
            step=step,
            key=slider_key
        )

        # Update session_state with slider values
        st.session_state[min_key] = selected_range[0]
        st.session_state[max_key] = selected_range[1]
        # Apply filter
        df_filtered = df[(df[selected_var] >= selected_range[0]) &
                         (df[selected_var] <= selected_range[1])]

    return selected_var, df_filtered

def render_thresholds_panel(df, label, thresholds_conditions):
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

    if miffi_labels:
        st.markdown("##### Miffi Labels Present in Dataset:")
        st.write(", ".join(map(str, miffi_labels)))
    else:
        st.markdown("##### Thresholds used for filtering:")
        for col, val in threshold_values.items():
            # st.metric(label=col, value=f"{val}")
            # Usar HTML para reducir el tamaño visual
            st.markdown(f"""
                <div style='text-align: center; font-size: 1em;'>
                    <em>{col}</em><br>
                    <span>{val}</span>
                </div>
            """, unsafe_allow_html=True)

    st.markdown("##### \n")
    filter_option = st.radio("Show data for:", ["All", "Accepted", "Rejected"])
    return filter_option

def render_aggrid_table(df, columns_tag, thresholds_conditions):
    df_filtered = df[[col for col in ALL_COLUMNS[columns_tag] if col in df.columns]].copy()
    visible_cols = [col for col in VISIBLE_COLUMNS[columns_tag] if col in df.columns]
    hidden_cols = [col for col in df_filtered.columns if col not in visible_cols]

    gb = GridOptionsBuilder.from_dataframe(df_filtered)
    gb.configure_selection('single')

    gb.configure_column(
        "micName",
        headerName="Micrograph",
        valueFormatter=JsCode("""
            function(params) {
                if (!params.value) return '';
                return params.value.split('/').slice(-1)[0];
            }
        """)
    )

    configure_aggrid_threshold_styles(gb, thresholds_conditions)

    for col in hidden_cols:
        gb.configure_column(col, hide=True)

    grid_response = AgGrid(
        df_filtered,
        gridOptions=gb.build(),
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        allow_unsafe_jscode=True,
        enable_enterprise_modules=False,
        height=400,
        theme="streamlit",
    )

    selected_rows = pd.DataFrame(grid_response['selected_rows'])

    if not selected_rows.empty:
        mic_path = selected_rows.iloc[0]["micName"]
        st.success(f"Selected: {mic_path}")
        return mic_path

    return ""

def render_micrograph_viewer(mic_path: str):
    st.markdown("##### Micrograph Viewer")

    if mic_path and os.path.exists(mic_path):
        try:
            with mrcfile.open(mic_path, permissive=True) as mrc:
                data = mrc.data

            image = (data - np.min(data)) / (np.max(data) - np.min(data))
            st.image(image, caption=os.path.basename(mic_path), clamp=True)
        except Exception as e:
            st.error(f"Error loading micrograph: {e}")
    else:
        st.info("No micrograph selected or file not found.")

def render_filter_view(df: pd.DataFrame, label: str, title: str, columns_tag: str, thresholds_conditions: dict,
                       plot_callback: Callable[[pd.DataFrame], str]):

    container_top_left, micrograph_viewer_panel = st.columns(2)
    with container_top_left:
        panel_distribution, panel_thresholds_button = st.columns(2)

    tab_table, tab_plots = st.tabs(["Dataframe", "Plots"])

    if label not in df.columns:
        st.warning(f'No views for the {title} yet')
        return

    with panel_distribution:
        plot_filters_distribution(df, label, title)

    with panel_thresholds_button:
        render_one_filter_statistics(df, columns_tag)
        filter_option = render_thresholds_panel(df, label, thresholds_conditions)

    df_filtered = apply_filter_selection(df, label, filter_option)

    mic_path_table = ""
    mic_path_plot = ""

    with tab_table:
        mic_path_table = render_aggrid_table(df_filtered, columns_tag, thresholds_conditions)

    with tab_plots:
        mic_path_plot = plot_callback(df_filtered)

    mic_path = mic_path_plot or mic_path_table

    with micrograph_viewer_panel:
        render_micrograph_viewer(mic_path)

def render_scores_view(df):
    container_top_left, container_top_right = st.columns(2)

    with container_top_left:
        plot_discrepancy_matrix(df)

    with container_top_right:
        plot_variable_correlation_matrix(df)

    tab_scatter, = st.tabs(["Correlation tab"])

    # Numeric variables
    visible_vars_raw = VISIBLE_COLUMNS.values()
    flattened_unique = list(dict.fromkeys(chain.from_iterable(visible_vars_raw)))
    df = df[flattened_unique]
    exclude_columns = ['micName', 'movieName']
    numeric_df = df.drop(columns=[col for col in df.columns if col.startswith('Filter') or col in exclude_columns])
    numeric_columns = numeric_df.columns.tolist()

    with tab_scatter:
        plot_correlation_scatter(df, numeric_columns)

def filter_visible_vars(df, visible_vars):
    excluded_prefixes = ['Filter']
    excluded_columns = ['micName', 'movieId', 'timestamp']

    filtered_vars = [col for col in visible_vars if col not in excluded_columns
                     and not any(col.startswith(prefix) for prefix in excluded_prefixes)
                     and pd.api.types.is_numeric_dtype(df[col])]

    return filtered_vars

# ------------------------------------- GENERAL PLOTS --------------------------------------
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

def plot_interactive_histogram(df, selected_var, default_bins=20):
    nbinsx = st.sidebar.slider("Bins number:", min_value=2, max_value=50, value=default_bins, step=5)

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=df[selected_var],
        nbinsx=nbinsx,
        name=selected_var,
        hovertemplate='Count: %{y}<br>' + selected_var + ': %{x:.2f}<extra></extra>'
    ))

    fig.update_layout(
        title=f'Histogram of {selected_var}',
        xaxis_title=selected_var,
        yaxis_title='Frecuency',
        bargap=0.2
    )

    return fig

def plot_interactive_scatter(df, selected_var, thresholds_conditions, time_col='movieId'):
    fig = go.Figure()

    lower_threshold = None
    upper_threshold = None
    added_keys = set()  # To avoid duplicates

    if selected_var in thresholds_conditions:
        condition = thresholds_conditions[selected_var]
        matches = re.findall(r'(-)?(threshold\w+)', condition)

        for sign, threshold_key in matches:
            key_id = f"{sign or ''}{threshold_key}"
            if key_id in added_keys:
                continue  # To avoid duplicates

            if threshold_key in df.columns:
                threshold_values = df[threshold_key]
                y_values = -threshold_values if sign == '-' else threshold_values
                label = f"-{threshold_key}" if sign == '-' else threshold_key

                fig.add_trace(go.Scatter(
                    x=df[time_col],
                    y=y_values,
                    mode='lines',
                    line=dict(color='red', dash='dash'),
                    name=f"Threshold: {label}",
                    hoverinfo='skip'  # Hide hover for thresholds
                ))

                if sign == '-':
                    lower_threshold = y_values
                else:
                    upper_threshold = y_values

                added_keys.add(key_id)

    y_values = df[selected_var]
    out_of_bounds = pd.Series([False] * len(df))

    if lower_threshold is not None and upper_threshold is not None:
        out_of_bounds = (y_values < lower_threshold) | (y_values > upper_threshold)
    elif lower_threshold is not None:
        out_of_bounds = y_values < lower_threshold
    elif upper_threshold is not None:
        out_of_bounds = y_values > upper_threshold
    else:
        out_of_bounds = pd.Series([False] * len(df))

    # Evaluate conditions dynamically
    if selected_var in thresholds_conditions:
        condition = thresholds_conditions[selected_var]
        for index, value in y_values.items():
            eval_condition = condition.replace('value', str(value))
            eval_condition = re.sub(
                r'threshold(\w+)',
                lambda m: f'df.loc[{index}, "threshold{m.group(1)}"]',
                eval_condition
            )
            try:
                out_of_bounds[index] = not eval(eval_condition)
            except Exception:
                out_of_bounds[index] = True

    colors = ['red' if out else 'steelblue' for out in out_of_bounds]

    fig.add_trace(go.Scatter(
        x=df[time_col],
        y=y_values,
        mode='markers',
        marker=dict(color=colors, size=7),
        name=selected_var,
        hovertemplate=f"{time_col}: %{{x}}<br>{selected_var}: %{{y}}<extra></extra>"
    ))

    fig.update_layout(
        title=f"{selected_var} vs {time_col}",
        # xaxis_title=time_col,
        yaxis_title=selected_var,
        height=500,
        width=1300,
        showlegend=False
    )

    return fig

# ------------------------------------- ANALYSIS FUNCTIONS --------------------------------------
def calculate_discrepancy_matrix(df, metric='discrepancy'):
    filter_columns = [col for col in df.columns if col.startswith('Filter')]
    matrix = np.zeros((len(filter_columns), len(filter_columns)))

    for i, col1 in enumerate(filter_columns):
        for j, col2 in enumerate(filter_columns):
            if i != j:
                if metric == 'discrepancy':
                    discrepancy = np.sum(df[col1] != df[col2])
                    matrix[i, j] = discrepancy / len(df) * 100
                elif metric == 'jaccard':
                    intersection = np.sum(df[col1] & df[col2])
                    union = np.sum(df[col1] | df[col2])
                    matrix[i, j] = intersection / union * 100
                elif metric == 'correlation':
                    col1_array = df[col1].fillna(False).astype(int).to_numpy()
                    col2_array = df[col2].fillna(False).astype(int).to_numpy()
                    correlation = np.corrcoef(col1_array, col2_array)[0, 1]
                    matrix[i, j] = correlation * 100
            else:
                if metric == 'discrepancy':
                    matrix[i, j] = 0
                elif metric == 'jaccard':
                    matrix[i, j] = 100
                elif metric == 'correlation':
                    matrix[i, j] = 100

    return filter_columns, matrix

def plot_matrix(filter_columns, matrix, metric):
    fig, ax = plt.subplots(figsize=(8, 6))
    if metric == 'discrepancy':
        sns.heatmap(matrix, xticklabels=filter_columns, yticklabels=filter_columns,
                    cmap='RdYlGn_r', annot=True, fmt=".2f", ax=ax, vmin=0, vmax=100)
    elif metric == 'jaccard':
        sns.heatmap(matrix, xticklabels=filter_columns, yticklabels=filter_columns,
                    cmap='RdYlGn', annot=True, fmt=".2f", ax=ax, vmin=0, vmax=100)
    elif metric == 'correlation':
        sns.heatmap(matrix, xticklabels=filter_columns, yticklabels=filter_columns,
                    cmap='RdYlGn', annot=True, fmt=".2f", ax=ax, vmin=-100, vmax=100)

    ax.set_title(f'{metric.capitalize()} Matrix (%)')
    ax.set_xlabel('Filters')
    ax.set_ylabel('Filters')
    plt.tight_layout()
    return fig

def calculate_correlation_matrix_variables(df):
    # Exclude specific columns and those starting with 'Filter'
    visible_vars_raw = VISIBLE_COLUMNS.values()
    flattened_unique = list(dict.fromkeys(chain.from_iterable(visible_vars_raw)))
    df = df[flattened_unique]
    exclude_columns = ['movieId', 'micName', 'movieName', 'creationTime']
    numeric_df = df.drop(columns=[col for col in df.columns if col.startswith('Filter') or col in exclude_columns])

    # Drop columns with all NaN values
    numeric_df = numeric_df.dropna(axis=1, how='all')

    # Convert categorical 'miffiLabel' to binary variables
    if 'miffiLabel' in numeric_df.columns:
        numeric_df = pd.get_dummies(numeric_df, columns=['miffiLabel'], drop_first=True)

    # Calculate correlation matrix
    correlation_matrix = numeric_df.corr()

    return correlation_matrix, numeric_df.columns

def plot_correlation_matrix_variables(correlation_matrix, columns, show_legend=False):
    fig, ax = plt.subplots(figsize=(14, 10))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, ax=ax)

    # Color the labels based on their group
    group_colors = {}
    for group, variables in VISIBLE_COLUMNS.items():
        color = next(ax._get_lines.prop_cycler)['color']
        for var in variables:
            if var in columns:
                group_colors[var] = color

    # Apply colors to x and y tick labels
    for label in ax.get_xticklabels():
        label.set_color(group_colors.get(label.get_text(), 'black'))
    for label in ax.get_yticklabels():
        label.set_color(group_colors.get(label.get_text(), 'black'))

    # Optionally add a legend
    if show_legend:
        patches = [mpatches.Patch(color=color, label=group) for group, color in group_colors.items()]
        ax.legend(handles=patches, title="Filter Groups", bbox_to_anchor=(1.05, 1), loc='upper left')

    ax.set_title('Correlation Matrix')
    plt.tight_layout()
    return fig

def plot_scatter(df, x, y):
    fig = px.scatter(df, x=x, y=y, title=f'Scatter plot of {x} vs {y}')
    fig.update_layout(
        title=dict(font=dict(size=12)),
        xaxis=dict(title=dict(font=dict(size=10))),
        yaxis=dict(title=dict(font=dict(size=10))),
        margin=dict(l=20, r=20, t=30, b=20)
    )
    return fig

# ------------------------------------- FILTERS VIEWS --------------------------------------
def main_view(df):
    render_filter_statistics(df)
    render_general_dataframe(df)

def dose_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'diffDosePerAngstrom2': "-thresholdPercentageDiff < value < thresholdPercentageDiff",
    }
    visible_vars_raw = VISIBLE_COLUMNS[columns_tag]
    visible_vars = filter_visible_vars(df, visible_vars_raw)
    selected_var, df_filtered = filter_slider_manager(df, visible_vars)

    def plot_tab(df_current):
        view_option = st.selectbox("", ["Interactive Scatter", "Interactive Histogram", "Others"])
        mic_path = ""
        if view_option == "Interactive Scatter":
            # time_col = 'movieId' if df_current['creationTime'].nunique() <= 1 else 'creationTime' esto hace que parpadee
            fig_scatter = plot_interactive_scatter(df_current, selected_var, thresholds_conditions)
            click_data = plotly_events(fig_scatter, click_event=True, hover_event=False)
            if click_data:
                point_index = click_data[0]['pointIndex']
                if 0 <= point_index < len(df_current):
                    row = df_current.iloc[point_index]
                    mic_path = row.get('micName', "")
        elif view_option == "Interactive Histogram":
            fig_hist = plot_interactive_histogram(df_current, selected_var)
            st.plotly_chart(fig_hist, use_container_width=True)
        elif view_option == "Others":
            st.write("Different filters plots")

        return mic_path

    render_filter_view(df_filtered, label, title, columns_tag, thresholds_conditions, plot_tab)

def drift_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'maxMovieShift': "value < thresholdMaxMovieShift",
        'maxFrameShift': "value < thresholdMaxFrameShift",
    }
    visible_vars_raw = VISIBLE_COLUMNS[columns_tag]
    visible_vars = filter_visible_vars(df, visible_vars_raw)
    selected_var, df_filtered = filter_slider_manager(df, visible_vars)

    def plot_tab(df_current):
        view_option = st.selectbox("", ["Interactive Scatter", "Interactive Histogram", "Others"])
        mic_path = ""
        if view_option == "Interactive Scatter":
            # time_col = 'movieId' if df_current['creationTime'].nunique() <= 1 else 'creationTime' esto hace que parpadee
            fig_scatter = plot_interactive_scatter(df_current, selected_var, thresholds_conditions)
            click_data = plotly_events(fig_scatter, click_event=True, hover_event=False)
            if click_data:
                point_index = click_data[0]['pointIndex']
                if 0 <= point_index < len(df_current):
                    row = df_current.iloc[point_index]
                    mic_path = row.get('micName', "")
        elif view_option == "Interactive Histogram":
            fig_hist = plot_interactive_histogram(df_current, selected_var)
            st.plotly_chart(fig_hist, use_container_width=True)
        elif view_option == "Others":
            st.write("Different filters plots")

        return mic_path

    render_filter_view(df_filtered, label, title, columns_tag, thresholds_conditions, plot_tab)

def tilt_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'tiltMeanCorrelation': "thresholdMeanCorrelation < value",
        'tiltStdCorrelation': "value < thresholdStdCorrelation",
    }

    visible_vars_raw = VISIBLE_COLUMNS[columns_tag]
    visible_vars = filter_visible_vars(df, visible_vars_raw)
    selected_var, df_filtered = filter_slider_manager(df, visible_vars)

    def plot_tab(df_current):
        view_option = st.selectbox("", ["Interactive Scatter", "Interactive Histogram", "Others"])
        mic_path = ""

        if view_option == "Interactive Scatter":
            # time_col = 'movieId' if df_current['creationTime'].nunique() <= 1 else 'creationTime' esto hace que parpadee
            fig_scatter = plot_interactive_scatter(df_current, selected_var, thresholds_conditions)
            click_data = plotly_events(fig_scatter, click_event=True, hover_event=False)
            if click_data:
                point_index = click_data[0]['pointIndex']
                if 0 <= point_index < len(df_current):
                    row = df_current.iloc[point_index]
                    mic_path = row.get('micName', "")
        elif view_option == "Interactive Histogram":
            fig_hist = plot_interactive_histogram(df_current, selected_var)
            st.plotly_chart(fig_hist, use_container_width=True)
        elif view_option == "Others":
            st.write("Different filters plots")

        return mic_path

    render_filter_view(df_filtered, label, title, columns_tag, thresholds_conditions, plot_tab)

def miffi_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'miffiLabel': "Filter Miffi"
    }
    visible_vars_raw = VISIBLE_COLUMNS[columns_tag]
    selected_var = 'miffiLabel'

    def plot_tab(df_current):
        view_option = st.selectbox("", ["Interactive Scatter", "Interactive Histogram", "Others"])
        mic_path = ""
        if view_option == "Interactive Scatter":
            pass
        elif view_option == "Interactive Histogram":
            fig_hist = plot_interactive_histogram(df_current, selected_var)
            st.plotly_chart(fig_hist, use_container_width=True)
        elif view_option == "Others":
            st.write("Different filters plots")

        return mic_path

    render_filter_view(df, label, title, columns_tag, thresholds_conditions, plot_tab)

def ctf_view(df, label, title, columns_tag):
    thresholds_conditions = {
        'defocusU': "thresholdMinDefocus < value < thresholdMaxDefocus",
        'defocusV': "thresholdMinDefocus < value < thresholdMaxDefocus",
        'astigmatismPercentage': "value < thresholdAstigmatismPercentage",
        'resolution': "value < thresholdResolution",
        'consensusResolution': "value < thresholdConsensusResolution",
    }

    visible_vars_raw = VISIBLE_COLUMNS[columns_tag]
    visible_vars = filter_visible_vars(df, visible_vars_raw)
    selected_var, df_filtered = filter_slider_manager(df, visible_vars)

    def plot_tab(df_current):
        view_option = st.selectbox("", ["Interactive Scatter", "Interactive Histogram", "Others"])
        mic_path = ""
        if view_option == "Interactive Scatter":
            # time_col = 'movieId' if df_current['creationTime'].nunique() <= 1 else 'creationTime' esto hace que parpadee
            fig_scatter = plot_interactive_scatter(df_current, selected_var, thresholds_conditions)
            click_data = plotly_events(fig_scatter, click_event=True, hover_event=False)
            if click_data:
                point_index = click_data[0]['pointIndex']
                if 0 <= point_index < len(df_current):
                    row = df_current.iloc[point_index]
                    mic_path = row.get('micName', "")
        elif view_option == "Interactive Histogram":
            fig_hist = plot_interactive_histogram(df_current, selected_var)
            st.plotly_chart(fig_hist, use_container_width=True)
        elif view_option == "Others":
            st.write("Different filters plots")

        return mic_path

    render_filter_view(df_filtered, label, title, columns_tag, thresholds_conditions, plot_tab)

def scores_view(df):
    st.sidebar.title("Plotting options")
    render_scores_view(df)

# ----------------------- SCORES VIEWS PLOTS --------------------------
def plot_discrepancy_matrix(df):
    st.markdown("##### Filter discrepancy matrix")
    metric = st.selectbox("Select metric", options=["discrepancy", "jaccard", "correlation"])
    filter_columns, matrix = calculate_discrepancy_matrix(df, metric)
    fig = plot_matrix(filter_columns, matrix, metric)
    st.pyplot(fig)

def plot_variable_correlation_matrix(df):
    st.markdown("##### Variable correlation matrix")

    if "variable_corr_matrix" not in st.session_state:
        st.session_state.variable_corr_matrix = None
        st.session_state.variable_corr_columns = None

    if st.button("Calculate variable correlation matrix"):
        correlation_matrix, columns = calculate_correlation_matrix_variables(df)
        st.session_state.variable_corr_matrix = correlation_matrix
        st.session_state.variable_corr_columns = columns

    if st.session_state.variable_corr_matrix is not None:
        fig = plot_correlation_matrix_variables(
            st.session_state.variable_corr_matrix,
            st.session_state.variable_corr_columns,
            show_legend=False
        )
        st.pyplot(fig)

def plot_correlation_scatter(df, numeric_columns):
    st.markdown("##### Correlation Scatter")

    if "scatter_plot_triggered" not in st.session_state:
        st.session_state.scatter_plot_triggered = False
    if "scatter_x" not in st.session_state:
        st.session_state.scatter_x = numeric_columns[0]
    if "scatter_y" not in st.session_state:
        st.session_state.scatter_y = numeric_columns[1]

    st.session_state.scatter_x = st.sidebar.selectbox("Select X-axis", numeric_columns, index=numeric_columns.index(st.session_state.scatter_x))
    st.session_state.scatter_y = st.sidebar.selectbox("Select Y-axis", numeric_columns, index=numeric_columns.index(st.session_state.scatter_y))

    if st.sidebar.button("Plot"):
        st.session_state.scatter_plot_triggered = True

    if st.session_state.scatter_plot_triggered:
        fig = plot_scatter(df, st.session_state.scatter_x, st.session_state.scatter_y)
        st.plotly_chart(fig, use_container_width=True)

# ------------------------------------------- MAIN FUNCTIONS -------------------------------------------------
def setStyle():
    st.markdown("""
        <style>    
            /* Fondo general del dashboard */
            .stApp {
            background-color: white;
            }
            
            /* Fondo del sidebar */
            section[data-testid="stSidebar"] {
                background-color: #1f3b57; /* Azul oscuro */
            }

            /* Fondo de los widgets */
            section[data-testid="stSidebar"] .stSelectbox,
            section[data-testid="stSidebar"] .stSlider,
            section[data-testid="stSidebar"] .stNumberInput,
            section[data-testid="stSidebar"] .stCheckbox,
            section[data-testid="stSidebar"] .stTextInput {
                background-color: #e6edf5; /* Azul muy claro */
                border-radius: 8px;
                padding: 10px;
            }

            /* Texto dentro de los widgets */
            section[data-testid="stSidebar"] input,
            section[data-testid="stSidebar"] select,
            section[data-testid="stSidebar"] textarea,
            section[data-testid="stSidebar"] .css-1cpxqw2, /* texto dentro de sliders */
            section[data-testid="stSidebar"] .st-bx,        /* texto dentro de selectbox */
            section[data-testid="stSidebar"] .st-cz {
                color: #0a2540 !important; /* Azul oscuro */
            }

            /* Títulos en el sidebar (como los de st.markdown("### ...")) */
            section[data-testid="stSidebar"] h1,
            section[data-testid="stSidebar"] h2,
            section[data-testid="stSidebar"] h3,
            section[data-testid="stSidebar"] h4,
            section[data-testid="stSidebar"] h5,
            section[data-testid="stSidebar"] h6 {
                color: white; /* Blanco */
            }

            /* Etiquetas de inputs */
            section[data-testid="stSidebar"] label {
                color: #1f3b57; /* Azul oscuro */
            }
            
            /* Títulos principales */
            h1, h2, h3, h4, h5, h6 {
                color: #1f3b57; /* Azul oscuro */
            }
            
        </style>
    """, unsafe_allow_html=True)

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

def load_and_validate_data(path):
    if not path or not os.path.exists(path):
        st.error("The CSV metadata file does not exist.")
        return None

    df = load_data(path)
    if df.empty:
        st.warning("The CSV is empty or could not be loaded.")
        return None

    return df

VIEW_FUNCTIONS = {
    'Main View': main_view,
    'Dose': lambda df: dose_view(df, 'Filter DoseAnalysis', 'Dose Analysis Filter', 'dose'),
    'Drift': lambda df: drift_view(df, 'Filter MaxShift', 'Drift Analysis Filter', 'drift'),
    'Tilt': lambda df: tilt_view(df, 'Filter TiltAnalysis', 'Tilt Analysis Filter', 'tilt'),
    'Miffi': lambda df: miffi_view(df, 'Filter Miffi', 'Miffi Filter', 'miffi'),
    'CTF': lambda df: ctf_view(df, 'Filter CTFConsensus', 'CTF Filter', 'ctf'),
    'Scores View': lambda df: scores_view(df),
}

def render_view(view_name, df):
    view_func = VIEW_FUNCTIONS.get(view_name)
    if view_func:
        view_func(df)
    else:
        st.error("Unknown view selected.")

def main():
    st.set_page_config(page_title="Quality Monitor", layout="wide")
    setStyle()
    st.title("Live Quality Metrics Monitor")
    st_autorefresh(interval=60_000, limit=None, key="quality_monitor_refresh")  # Auto-refresh every 30 seconds (30,000 milliseconds)
    monitor_path = os.environ.get("STREAMLIT_MONITOR_PATH") # .csv metadata

    df = load_and_validate_data(monitor_path)
    if df is None:
        return

    # --------------  BUILD DASHBOARD ------------------------------
    # View selector ('Main View', 'Dose', 'Drift', 'Tilt', 'Miffi', 'CTF', 'Scores View')
    st.sidebar.markdown("## Object Views")
    selected_view = st.sidebar.selectbox('Select option:', VIEWS)

    # --- Entries range filter ---
    st.sidebar.markdown("### Entries range filter")
    enable_range_filter = st.sidebar.checkbox("Activate range filter")
    entry_range = (0, len(df))

    if enable_range_filter:
        sort_col = 'movieId' if df['creationTime'].nunique() <= 1 else 'creationTime'
        df = df.sort_values(by=sort_col).reset_index(drop=True)

        max_index = len(df) - 1
        default_start = max(0, max_index - 99)
        entry_range = st.sidebar.slider("Select entries range:", min_value = 0, max_value = max_index,
                                        value = (default_start, max_index), step = 1)
        df = df.iloc[entry_range[0]:entry_range[1] + 1]

    render_view(selected_view, df)


if __name__ == "__main__":
    main()
