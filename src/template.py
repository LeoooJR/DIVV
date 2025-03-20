from jinja2 import Environment, FileSystemLoader
from pandas import isna

class Report:

    __slots__ = ('vcfs', 'infos', 'cmd', 'view', 'plots', 'prefix', 'summary')

    def __init__(self, vcfs: list[str], prefix: str, cmd: str, infos: dict, view: object, plots: dict, summary: object = None):

        self.vcfs = vcfs

        self.cmd = cmd

        self.infos = infos
        
        self.view = view

        plots[vcfs[1]].dark()

        self.plots = {vcfs[0]:plots[vcfs[0]].as_html(),
                      vcfs[1]:plots[vcfs[1]].as_html(),
                      "common":plots["common"].as_html()}

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

        html = template.render(vcfs=self.vcfs, cmd=self.cmd, infos=self.infos, view=self.view, summary=self.summary, plots=self.plots)

        with open(f"{self.prefix}.html",'w') as f:
            f.writelines(html)