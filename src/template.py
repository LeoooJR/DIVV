from glob import glob
import os
from jinja2 import Environment, FileSystemLoader
from pandas import isna
from shutil import copy, copytree

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
        
        ressources = os.path.join(os.path.dirname(os.path.abspath(__file__)),'templates/')
        
        env = Environment(loader=FileSystemLoader(ressources))

        env.filters['is_nan'] = is_nan

        template = env.get_template("template.html")

        html = template.render(vcfs=self.vcfs, cmd=self.cmd, infos=self.infos, view=self.view, summary=self.summary, plots=self.plots)

        with open(f"{self.prefix}.html",'w') as f:
            f.writelines(html)

        for f in glob(f"{ressources}*.css"):
            copy(f,os.getcwd())
        copytree(os.path.join(ressources,'statics'),os.path.join(os.getcwd(),'statics'), dirs_exist_ok=True)