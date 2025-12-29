import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import io

# Set page configuration
st.set_page_config(page_title="Metabolomics Heatmap Visualizer", layout="wide")

def main():
    st.title("🧪 Metabolomics Heatmap Visualizer")
    st.markdown("""
    Upload your metabolomics CSV file, map the columns, and customize your visualization.
    This app supports log transformation, normalization, and **vector graphics export (SVG/PDF)**.
    """)

    # Sidebar for Configuration
    st.sidebar.header("1. Data Upload")
    uploaded_file = st.sidebar.file_uploader("Upload your CSV file", type=["csv"])

    if uploaded_file is not None:
        try:
            # Load data
            df = pd.read_csv(uploaded_file)
            st.sidebar.success("File uploaded successfully!")

            # 2. Column Selection
            st.sidebar.header("2. Column Mapping")
            all_columns = df.columns.tolist()

            # Default selections based on your specific file structure
            feature_id_col = st.sidebar.selectbox("Select Feature ID column", all_columns, index=1 if len(all_columns) > 1 else 0)
            compound_name_col = st.sidebar.selectbox("Select Compound Name column", all_columns, index=min(36, len(all_columns)-1))
            
            # Suggesting C-T columns if they exist as a default range (indices 2 to 20)
            default_samples = all_columns[2:20] if len(all_columns) >= 20 else []
            sample_cols = st.sidebar.multiselect(
                "Select Sample columns (Intensities)", 
                all_columns, 
                default=default_samples
            )

            if not sample_cols:
                st.warning("Please select at least one sample column in the sidebar.")
            else:
                # 3. Processing Options
                st.sidebar.header("3. Processing Options")
                use_log10 = st.sidebar.checkbox("Apply Log10 Transformation", value=True)
                use_normalization = st.sidebar.checkbox("Apply Normalization (Center Scaling)", value=False, 
                                                   help="Subtracts mean and divides by standard deviation for each compound.")
                
                # 4. Visualization Options
                st.sidebar.header("4. Heatmap Settings")
                color_theme = st.sidebar.selectbox(
                    "Color Palette",
                    ["YlGnBu", "Viridis", "Plasma", "Inferno", "Magma", "RdBu_r", "Rocket", "Mako"]
                )
                
                font_size = st.sidebar.slider("Label Font Size", 5, 20, 10)

                # Data Processing
                plot_df = df[[compound_name_col, feature_id_col] + sample_cols].copy()
                
                # Handle duplicates/NAs
                plot_df[compound_name_col] = plot_df[compound_name_col].fillna("Unknown")
                plot_df['Display_Name'] = plot_df[compound_name_col].astype(str) + " (" + plot_df[feature_id_col].astype(str) + ")"
                plot_df = plot_df.set_index('Display_Name')
                plot_df = plot_df[sample_cols]

                # Convert to numeric
                plot_df = plot_df.apply(pd.to_numeric, errors='coerce').fillna(0)

                # Apply Log10
                if use_log10:
                    plot_df = np.log10(plot_df + 1)

                # Apply Normalization
                if use_normalization:
                    row_mean = plot_df.mean(axis=1)
                    row_std = plot_df.std(axis=1).replace(0, 1)
                    plot_df = plot_df.sub(row_mean, axis=0).div(row_std, axis=0)

                # Create Heatmap
                fig = px.imshow(
                    plot_df,
                    labels=dict(x="Samples", y="Compounds", color="Value"),
                    x=plot_df.columns,
                    y=plot_df.index,
                    color_continuous_scale=color_theme,
                    aspect="auto",
                    title=f"Heatmap: {'Normalized ' if use_normalization else ''}{'Log10 ' if use_log10 else ''}Intensities"
                )

                fig.update_layout(
                    height=max(600, len(plot_df) * 15 + 200),
                    xaxis_nticks=len(sample_cols),
                    yaxis_autorange='reversed',
                    font=dict(size=font_size)
                )

                # --- DISPLAY & EXPORT ---
                st.subheader("📊 Heatmap Visualization")
                st.plotly_chart(fig, use_container_width=True)

                st.divider()
                st.subheader("💾 Export Heatmap (Vector Graphics)")
                exp_col1, exp_col2, exp_col3 = st.columns(3)

                with st.spinner("Preparing export files..."):
                    try:
                        svg_data = fig.to_image(format="svg", engine="kaleido")
                        pdf_data = fig.to_image(format="pdf", engine="kaleido")
                        png_data = fig.to_image(format="png", width=1200, height=max(800, len(plot_df)*20), scale=2)

                        with exp_col1:
                            st.download_button(label="Download SVG (Vector)", data=svg_data, file_name="heatmap.svg", mime="image/svg+xml")
                        with exp_col2:
                            st.download_button(label="Download PDF (Vector)", data=pdf_data, file_name="heatmap.pdf", mime="application/pdf")
                        with exp_col3:
                            st.download_button(label="Download PNG (High-Res)", data=png_data, file_name="heatmap.png", mime="image/png")
                    except Exception as e:
                        st.error(f"Export Error: {e}")

                st.divider()
                st.subheader("📄 Processed Data Preview")
                st.dataframe(plot_df.head(10))
                csv = plot_df.to_csv().encode('utf-8')
                st.download_button("Download Processed CSV", csv, "processed_metabolomics_data.csv", "text/csv")

        except Exception as e:
            st.error(f"Critical Error: {e}")
    else:
        st.info("Welcome! Please upload your metabolomics CSV file in the sidebar to begin.")

    # Sidebar Footer (Logo and Website)
    st.sidebar.markdown("---")
    # Using the raw URL for the logo so it renders properly in Streamlit
    logo_url = "https://raw.githubusercontent.com/trikaloudis/aquomixlab_heatmaps/main/Aquomixlab%20Logo%20v2%20white%20font.jpg"
    st.sidebar.image(logo_url, use_container_width=True)
    st.sidebar.markdown(
        """
        <div style="text-align: center; margin-top: -10px;">
            <a href="https://www.aquomixlab.com/" target="_blank" style="text-decoration: none; color: #4A90E2; font-weight: bold; font-family: sans-serif;">
                Visit Aquomixlab Website
            </a>
        </div>
        """, 
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
