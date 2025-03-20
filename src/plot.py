import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from memory_profiler import profile
import numpy
import plotly.express as px
import plotly.graph_objs as go
import plotly.io as pio
from pandas import Series, DataFrame, concat
import pprint
import utils

class Plot:

    __slots__ = ('fig')

    def __init__(self, fig: object):
        
        self.fig = fig

    def __str__(self):
        pass

    @staticmethod
    def set_layout(fig: object):

        fig.update_layout(
            font=dict(
                size = 15
            ),
            title=dict(
                x = 0.5,
                xanchor= 'center',
                yanchor= 'top'
            ),
            barmode='relative',
            hovermode="x unified",
            legend=dict(
                traceorder="normal",
                itemclick="toggle",
                itemdoubleclick="toggleothers"
            )
        )
    
class PlotLibrary:

    __slots__ = ('file', 'plots')

    def __init__(self, file: str = None):
        
        self.file = file
        self.plots = []

    def __str__(self):
        return f"Library with {len(self.plots)} plots for file {self.file}"

    def save(self, plot: Plot):
        
        self.plots.append(plot)

    def as_html(self):
         return list(map(lambda p: pio.to_html(p.fig, full_html=False, include_plotlyjs=False),self.plots))

    def barplot(self, data, nominal: str, y: str, color: str, title: str, prefix: str) -> Plot:

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
                    title=title)
            
        fig.update_xaxes(ticklabelstep=1)

        plot = Plot(fig=fig)

        Plot.set_layout(plot.fig)

        self.save(plot)

    def histogram(self, df: DataFrame, x: str, color: str, title: str, prefix: str) -> Plot:

        df.dropna(axis=0, how='any', inplace=True)

        fig = px.histogram(data_frame=df, 
                           x=x, 
                           title=title)

        plot = Plot(fig=fig)

        Plot.set_layout(plot.fig)

        self.save(plot)

    def boxplot(self, df: DataFrame, nominal: str, y: str, color: str, title: str, prefix: str) -> Plot:

        fig = px.box(data_frame=df,
                        x=nominal,
                        y=y,
                        color=color,
                        title=title)
            
        plot = Plot(fig=fig)

        Plot.set_layout(plot.fig)
        
        self.save(plot)

    def venn(self, sizes: tuple[int], labels: list[str] = None) -> Plot:

        NSETS = 2
        NSUBSETS = 3
        PADDING = 0.2
        
        v = venn2(sizes,labels)

        plt.close()

        colors: list = ['Tomato', 'LightGreen']

        shapes: list = [go.layout.Shape(type="circle",
                         xref="x",
                         yref="y", 
                         x0=v.centers[i].x - v.radii[i], 
                         y0=v.centers[i].y - v.radii[i], 
                         x1=v.centers[i].x + v.radii[i], 
                         y1=v.centers[i].y + v.radii[i], 
                         fillcolor=colors[i], 
                         line_color=colors[i], 
                         opacity=0.75) for i in range(0,NSETS)]
        
        annotations: list = [go.layout.Annotation(xref="x",
                                                  yref="y", 
                                                  x= v.set_labels[i].get_position()[0], 
                                                  y=v.set_labels[i].get_position()[1], 
                                                  text= v.set_labels[i].get_text(), 
                                                  showarrow=False) for i in range(0,NSETS)]
        
        annotations.extend(go.layout.Annotation(xref="x",
                                                yref="y", 
                                                x= v.subset_labels[i].get_position()[0], 
                                                y=v.subset_labels[i].get_position()[1], 
                                                text= v.subset_labels[i].get_text(), 
                                                showarrow=False) for i in range(0,NSUBSETS))
        
        xmin = min(v.centers[i].x - v.radii[i] for i in range(0,NSETS)) - PADDING
        xmax = max(v.centers[i].x + v.radii[i] for i in range(0,NSETS)) + PADDING
        ymin = min(v.centers[i].y - v.radii[i] for i in range(0,NSETS)) - PADDING
        ymax = max(v.centers[i].y + v.radii[i] for i in range(0,NSETS)) + PADDING
        
        fig = go.Figure()

        fig.update_xaxes(range=[xmin,xmax], showticklabels=False, ticklen=0)

        fig.update_yaxes(range=[ymin,ymax], showticklabels=False, ticklen=0, scaleanchor="x", scaleratio=1)

        fig.update_layout(
            plot_bgcolor='white',
            margin = dict(b = 0, l = 10, pad = 0, r = 10, t = 40),
            width=800, 
            height=400,
            shapes=shapes, 
            annotations=annotations,
        )

        plot = Plot(fig=fig)

        self.save(plot)

    def dark(self):
        list(map(lambda plot: plot.fig.update_layout(template="plotly_dark"),self.plots))

def visualization(file: str, stats: object):

    library = PlotLibrary(file=file)

    data = []

    chromosomes = stats.keys()

    for k in chromosomes:

        data.append({"Chromosome": k, "Type": "Indel", "Count": stats[k]["variant"]["indel"]["deletion"] + stats[k]["variant"]["indel"]["insertion"]})
        data.append({"Chromosome": k, "Type": "SNP", "Count": stats[k]["variant"]["snp"]["transition"] + stats[k]["variant"]["snp"]["transversion"]})
        data.append({"Chromosome": k, "Type": "SV", "Count": stats[k]["variant"]["sv"]})

    library.barplot(data, "Chromosome","Count", "Type", "Variant by chromosome", "VariantByChromosome")

    data = DataFrame([{"Chromosome": k, "Genotype": genotype, "Count": stats[k][code]} 
                    for k in chromosomes 
                    for genotype, code in [("Homozygous", "hom"), ("Heterozygous", "het")]])

    library.barplot(data, "Chromosome", "Count", "Genotype", "Genotype by chromosome","GenotypeByChromosome")

    data = DataFrame([{"Chromosome": k, "SNP": snp, "Count": stats[k]["variant"]["snp"][snp]} for k in chromosomes for snp in ["transition","transversion"]])

    library.barplot(data, "Chromosome", "Count", "SNP", "SNP Type by Chromosome", "SNPTypeByChrom")

    chromosome = list(chromosomes)[0]

    if stats[chromosome]["GQ"].size:
        pass

    if stats[chromosome]["depth"].size:

        with utils.suppress_warnings():       
            data = DataFrame(list(map(lambda k: Series([k, numpy.mean(stats[k]["depth"])], index=["Chromosome", "Depth"]), chromosomes)))
        
        library.barplot(data, "Chromosome", "Depth", color="Chromosome", title="Mean depth by chromosome", prefix="DepthByChromosomeBarPlot")

        data = list(map(lambda k: DataFrame({"Chromosome": [k] * stats[k]["depth"].size,
                                 "Depth": stats[k]["depth"].flatten()}), chromosomes))

        df = concat(data, ignore_index=True).astype({"Chromosome": "category", "Depth": "int"})

        library.boxplot(df, "Chromosome", "Depth", "Chromosome", "Depth by chromosome", "DepthByChromosomeBoxPlot")

        library.histogram(df, "Depth", None, "Depth", "DepthHist")

        data.clear()

    if stats[chromosome]["quality"].size:
        pass

    return library