import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import plotly.express as px
import plotly.graph_objs as go
import plotly.io as pio
from pandas import DataFrame, concat

class Plot:
    """ 
    A class representing a plot 

    Args:
        fig: The figure object
    """

    __slots__ = ('fig')

    def __init__(self, fig: object):
        
        self.fig = fig

        self.set_layout()

    def set_layout(self):

        self.fig.update_layout(
            font=dict(
                size = 15,
                family = "ui-sans-serif, system-ui, sans-serif, Apple Color Emoji, Segoe UI Emoji, Segoe UI Symbol, Noto Color Emoji",
                weight = 400
            ),
            title=dict(
                x = 0.5,
                xanchor= 'center',
                yanchor= 'top'
            ),
            hovermode="x unified",
            legend=dict(
                traceorder="normal",
                itemclick="toggle",
                itemdoubleclick="toggleothers"
            )
        )
    
    def __str__(self):
        pass
    
class PlotLibrary:
    """ A class representing a library of plots
    
    Args:
        file: The path to the file associated with the plots
    """

    __slots__ = ('file', 'plots')

    def __init__(self, file: str = None):
        
        self.file = file
        self.plots = []

    def __str__(self):
        return f"Library with {len(self.plots)} plots about {self.file}"

    def save(self, plot: Plot):
        """ Save the plot to the library """
        
        self.plots.append(plot)

    def as_html(self):
        """ Return the plots as HTML """

        return list(map(lambda p: pio.to_html(p.fig, full_html=False, include_plotlyjs=False, config={
            "displayModeBar":"hover",
            "displaylogo":False}),self.plots))

    def barplot(self, data, nominal: str, y: str, color: str, title: str, prefix: str) -> Plot:
        """ Create a bar plot """

        if isinstance(data, list):
            if isinstance(data[0],DataFrame):
                df = concat(data, ignore_index=True)
            elif isinstance(data[0], dict):
                df = DataFrame(data)
        else:
            df = data

        fig = px.bar(data_frame=df[df[y] > 0], 
                    x=nominal, 
                    y=y, 
                    color=color,
                    title=title,
                    barmode="relative")
        
        fig.update_traces(offsetgroup=None)
        
        # Set x-axis type to category for proper stacking
        fig.update_layout(xaxis_type='category')

        plot = Plot(fig=fig)

        self.save(plot)

    def histogram(self, df: DataFrame, x: str, color: str, title: str, prefix: str) -> Plot:
        """ Create a histogram plot """

        df.dropna(axis=0, how='any', inplace=True)

        fig = px.histogram(data_frame=df, 
                           x=x, 
                           title=title)
        
        fig.update_layout(yaxis=dict(
            title=dict(
                text="Count"
            )
        ))

        plot = Plot(fig=fig)

        self.save(plot)

    def boxplot(self, df: DataFrame, nominal: str, y: str, color: str, title: str, prefix: str) -> Plot:
        """ Create a box plot """
        fig = px.box(data_frame=df,
                        x=nominal,
                        y=y,
                        color=color,
                        title=title)
            
        plot = Plot(fig=fig)
        
        self.save(plot)

    def venn(self, sizes: tuple[int], labels: list[str] = None) -> Plot:
        """ Create a Venn diagram """

        NSETS: int = 2
        NSUBSETS: int = 3
        PADDING: float = 0.2
        
        v = venn2(sizes, labels)

        plt.close()

        colors: list[str] = ['#2C7A7B', '#822E4A']

        shapes: list[go.layout.Shape] = [go.layout.Shape(type="circle",
                         xref="x",
                         yref="y", 
                         x0=v.centers[i].x - v.radii[i], 
                         y0=v.centers[i].y - v.radii[i], 
                         x1=v.centers[i].x + v.radii[i], 
                         y1=v.centers[i].y + v.radii[i], 
                         fillcolor=colors[i], 
                         line_color=colors[i], 
                         opacity=0.75, 
                         name=labels[i]) for i in range(0,NSETS)]
                
        annotations: list[go.layout.Annotation] = [go.layout.Annotation(xref="x",
                                                  yref="y", 
                                                  x=v.set_labels[i].get_position()[0] - v.radii[i] if not i else v.set_labels[i].get_position()[0] + v.radii[i], 
                                                  y=v.set_labels[i].get_position()[1], 
                                                  text=v.set_labels[i].get_text(), 
                                                  showarrow=False) for i in range(0,NSETS)]
        
        annotations.extend(go.layout.Annotation(xref="x",
                                                yref="y", 
                                                x=v.subset_labels[i].get_position()[0], 
                                                y=v.subset_labels[i].get_position()[1], 
                                                text=v.subset_labels[i].get_text(), 
                                                showarrow=False) for i in range(0,NSUBSETS))
        
        xmin: float = min(v.centers[i].x - v.radii[i] for i in range(0,NSETS)) - PADDING
        xmax: float = max(v.centers[i].x + v.radii[i] for i in range(0,NSETS)) + PADDING
        ymin: float = min(v.centers[i].y - v.radii[i] for i in range(0,NSETS)) - PADDING
        ymax: float = max(v.centers[i].y + v.radii[i] for i in range(0,NSETS)) + PADDING
        
        fig: go.Figure = go.Figure()

        fig.update_xaxes(range=[xmin,xmax], showticklabels=False, ticklen=0)

        fig.update_yaxes(range=[ymin,ymax], showticklabels=False, ticklen=0, scaleanchor="x", scaleratio=1)

        fig.update_layout(
            paper_bgcolor='#fafafa',
            plot_bgcolor='#fafafa',
            margin = dict(b = 0, l = 10, pad = 0, r = 10, t = 40),
            width=800, 
            height=400,
            shapes=shapes, 
            annotations=annotations,
            title="Venn Diagram"
        )

        plot: Plot = Plot(fig=fig)

        self.save(plot)

    def dark(self) -> None:
        """ Set the plots to a dark theme """
        list(map(lambda plot: plot.fig.update_layout(template="plotly_dark"),self.plots))

    def light(self) -> None:
        """ Set the plots to a light theme """
        list(map(lambda plot: plot.fig.update_layout(template="plotly_light"),self.plots))