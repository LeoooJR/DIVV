from jinja2 import Environment, FileSystemLoader
from pandas import isna

class Report:

    __slots__ = ('vcfs', 'infos', 'data', 'plots', 'prefix', 'summary')

    def __init__(self, vcfs: list[str], prefix: str, infos: dict, df: object, plots: dict, summary: object = None):

        self.vcfs = vcfs

        self.infos = infos
        
        self.data = df

        self.plots = plots

        self.prefix = prefix

        self.summary = summary

    def __str__(self):
        pass

    def create(self):

        def is_nan(value):
            return isna(value)
        
        env = Environment(loader=FileSystemLoader('src/templates'))

        env.filters['is_nan'] = is_nan

        template = env.get_template("template.html")

        self.plots[self.vcfs[1]].dark()

        html = template.render(vcfs=self.vcfs, infos=self.infos, df=self.data, summary=self.summary, plots={self.vcfs[0]:self.plots[self.vcfs[0]].as_html(),
                                                                                                            self.vcfs[1]:self.plots[self.vcfs[1]].as_html()})

        with open(f"{self.prefix}.html",'w') as f:
            f.writelines(html)