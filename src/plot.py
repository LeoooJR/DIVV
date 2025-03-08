import numpy
import plotly.express as px
import pandas as pd
import pprint

class Plot:

    def __init__(self, type: str):
        pass

    def __str__(self):
        pass

def set_layout(fig: object) -> object:

    fig.update_layout(
        font=dict(
            size = 15
        ),
        title=dict(
            x = 0.5,
            xanchor= 'center',
            yanchor= 'top'
        )
    )

    return fig

def barplot(data: list, nominal: str, y: str, color: str, title: str, prefix: str):

    if isinstance(data[0],pd.DataFrame):
        df = pd.concat(data, ignore_index=True)
    elif isinstance(data[0],dict):
        df = pd.DataFrame(data)

    fig = px.bar(data_frame=df, 
                 x=nominal, 
                 y=y, 
                 color=color,
                 title=title)
    
    fig.update_xaxes(ticklabelstep=1)
    
    fig.update_layout(
        updatemenus=[
        dict(
            buttons=list([
                dict(
                    args=["type", "bar"],
                    label="Bar Chart",
                    method="restyle"
                )
            ]),
            direction="down",
        ),
    ]
    )

    set_layout(fig).write_html(f"{prefix}.html")

def histogram(data: list, x: str, color: str, title: str, prefix: str):

    if isinstance(data[0],pd.DataFrame):
        df = pd.concat(data, ignore_index=True)
    elif isinstance(data[0],dict):
        df = pd.DataFrame(data)

    fig = px.histogram(data_frame=df, x=x, title=title)

    set_layout(fig).write_html(f"{prefix}.html")

def boxplot(data: list, nominal: str, y: str, color: str, title: str, prefix: str):

    if isinstance(data[0],pd.DataFrame):
        df = pd.concat(data, ignore_index=True)
    elif isinstance(data[0],dict):
        df = pd.DataFrame(data)

    fig = px.box(data_frame=df,
                     x=nominal,
                     y=y,
                     color=color,
                     title=title)
    
    set_layout(fig).write_html(f"{prefix}.html")

def visualization(file: str, stats: object):

    data = []

    chromsomes = stats.keys()

    for k in chromsomes:

        data.append({"Chromosome": k, "Type": "Indel", "Count": stats[k]["variant"]["indel"]["deletion"] + stats[k]["variant"]["indel"]["insertion"]})
        data.append({"Chromosome": k, "Type": "SNP", "Count": stats[k]["variant"]["snp"]["transition"] + stats[k]["variant"]["snp"]["transversion"]})
        data.append({"Chromosome": k, "Type": "SV", "Count": stats[k]["variant"]["sv"]})

    barplot(data,"Chromosome","Count","Type", "Variant by chromosome", "VariantByChromosome")

    data.clear()

    for k in chromsomes:
        data.append({"Chromosome": k, "Genotype": "Homozygous", "Count": stats[k]["hom"]})
        data.append({"Chromosome": k, "Genotype": "Heterozygous", "Count": stats[k]["het"]})

    barplot(data,"Chromosome","Count","Genotype", "Genotype by chromosome","GenotypeByChromosome")

    data.clear()

    chromosome = list(chromsomes)[0]

    if stats[chromosome]["GQ"].size:
        pass

    if stats[chromosome]["depth"].size:

        for k in chromsomes:
            data.append({"Chromosome": k, "Depth": numpy.mean(stats[k]["depth"])})
        
        barplot(data, "Chromosome", "Depth", color="Chromosome", title="Depth by chromosome", prefix="DepthByChromosome")

        data.clear()

        for k in chromsomes:

            tmp = pd.DataFrame({"Chromosome": [k] * stats[k]["depth"].size,
                                "Depth": stats[k]["depth"].flatten()})

            data.append(tmp)

        boxplot(data, "Chromosome", "Depth", "Chromosome", "Depth by chromosome", "DepthByChromosome")

        histogram(data, "Depth", None, "Depth", "DepthHist")

        data.clear()

    if stats[chromosome]["quality"].size:
        pass