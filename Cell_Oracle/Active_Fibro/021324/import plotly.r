import plotly.express as px
import plotly.graph_objects as go

def plot_bubble_plot_2(df, title, output_file):

    # Create the bubble plot
    fig = px.scatter(
        df,
        x="TF",
        y="Status",
        size="Fraction",
        color="Median Differences",
        color_continuous_scale="Viridis",
        title=title,
        labels={
            "Fraction": "Fraction of Cells",
            "Median Differences": "Median Differences"
        },
        size_max=25,
    )

    # Add x-axis and y-axis lines and expand plot dimensions
    fig.update_layout(
        xaxis_title="TF",
        yaxis_title="Status",
        template="plotly_white",
        yaxis=dict(dtick=1),
        shapes=[
            # x-axis line
            dict(
                type="line",
                x0=-0.5,
                y0=-0.5,
                x1=len(df["TF"].unique()) - 0.5,
                y1=-0.5,
                line=dict(color="black", width=2)
            ),
            # y-axis line
            dict(
                type="line",
                x0=-0.5,
                y0=-0.5,
                x1=-0.5,
                y1=len(df["Status"].unique()) - 0.5,
                line=dict(color="black", width=2)
            )
        ],
    )

    # Shrink the colorbar legend
    fig.update_traces(
        marker=dict(
            colorbar=dict(
                thickness=15,
                len=0.5,
                title=dict(
                    text="Median Differences",
                    side="right"
                )
            )
        )
    )

    # Add a manual size legend (custom annotations or scatter points)
    max_size = 25  # Maximum bubble size (set earlier)
    size_legend_values = [0.1, 0.3, 0.6, 1.0]  # Example values to show as legends
    size_legend_markers = [v * max_size for v in size_legend_values]

    # Add scatter traces for the size legend
    for i, size in enumerate(size_legend_markers):
        fig.add_trace(
            go.Scatter(
                x=[-1.2],  # Position slightly outside the x-axis
                y=[len(df["Status"].unique()) - 0.5 - i],
                mode="markers+text",
                marker=dict(size=size, color="gray", opacity=0.6),
                text=[f"{size_legend_values[i]:.1f}"],  # Show size value
                textposition="middle right",
                showlegend=False
            )
        )

    # Add legend label
    fig.add_annotation(
        x=-1.4,
        y=len(df["Status"].unique()) - 0.2,
        text="Fraction of Cells (Bubble Size)",
        showarrow=False,
        font=dict(size=12)
    )

    # Write the figure to a file and show
    fig.write_image(output_file, format='pdf')
    fig.show()
