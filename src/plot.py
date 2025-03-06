import numpy
import plotly.express as px
import pandas as pd
import pprint

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

    fig.write_html(f"{prefix}.html")

def histogram(data: list, nominal: str, y: str, color: str, title: str, prefix: str):
    pass

def boxplot(data: list, nominal: str, y: str, color: str, title: str, prefix: str):
    pass

def visualization(stats: object):

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

    key = list(stats.keys())[0]

    if stats[key]["GQ"].size:
        pass

    if stats[key]["depth"].size:

        for k in chromsomes:
            data.append({"Chromosome": k, "Depth": numpy.mean(stats[k]["depth"])})
        
        barplot(data, "Chromosome", "Depth", color="Chromosome", title="Depth by chromosome", prefix="DepthByChromosome")

        data.clear()

        for k in chromsomes:

            tmp = pd.DataFrame({"Chromosome": [k] * stats[k]["depth"].size,
                                "Depth": stats[k]["depth"].flatten()})

            data.append(tmp)
        
        df = pd.concat(data, ignore_index=True)

        fig = px.box(data_frame=df,
                     x="Chromosome",
                     y="Depth")
        
        #fig.show()

        fig = px.histogram(data_frame=df, x="Depth")

        fig.show()

    if stats[key]["quality"].size:
        pass