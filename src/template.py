from glob import glob
import os
from jinja2 import Environment, FileSystemLoader
from pandas import isna, DataFrame
from shutil import copy, copytree
import zipfile

class Report:
    """
    A class to represent a report
        vcfs: The list of VCF files
        prefix: The prefix of the report
        cmd: The command used to generate the report
        infos: The information about the VCF files
        view: The view object as a DataFrame
        plots: The plots
        table: The benchmark table (truth metrics) object as a DataFrame
    """
    # Make use of __slots__ to avoid the creation of __dict__ and __weakref__ for each instance, reducing the memory footprint
    __slots__ = ('vcfs', 'infos', 'cmd', 'view', 'plots', 'prefix', 'table', 'out')

    def __init__(self, vcfs: list[str], prefix: str, cmd: str, infos: dict, view: dict, plots: dict, table: DataFrame = None, out: str = "archive"):

        # The list of VCF files
        self.vcfs = vcfs

        # The command used to generate the report
        self.cmd = cmd

        # The information about the VCF files
        self.infos = infos
        
        # The view object as a DataFrame
        self.view = view

        # Set one library of plots to dark mode
        plots[vcfs[1]].dark()

        # The plots as html
        self.plots = {vcfs[0]:plots[vcfs[0]].as_html(),
                      vcfs[1]:plots[vcfs[1]].as_html(),
                      "common":plots["common"].as_html()}

        # The prefix of the report
        self.prefix = prefix

        # The benchmark table (truth metrics) object as a DataFrame
        self.table = table

        # How to save the outputs
        self.out = out

    def __str__(self):
        pass

    def create(self):

        # Custom filter to check if a value is NaN with Jinja2
        def is_nan(value):
            return isna(value)
        
        # Path to look for the template
        ressources = os.path.join(os.path.dirname(os.path.abspath(__file__)),'templates')

        # Path to look for the CSS stylesheets
        stylesheets = os.path.join(ressources, "*.css")

        # Path to look for the statics
        statics = os.path.join(ressources, 'statics', "*.png")

        # HTML output
        output = f"{self.prefix}.html"
        
        # Create the environment
        env = Environment(loader=FileSystemLoader(ressources))

        # Add the custom filter to the environment
        env.filters['is_nan'] = is_nan

        # Load the template
        template = env.get_template("template.html")

        # Render the template
        html = template.render(vcfs=self.vcfs, cmd=self.cmd, infos=self.infos, view=self.view, table=self.table, plots=self.plots)

        if self.out == "dir":

            # Open data stream
            try:
                os.mkdir(f"{self.prefix}")
            except (FileNotFoundError, NotImplementedError, PermissionError):
                raise IOError("Report cannot be created.")

            with open(os.path.join(self.prefix, output),'w') as f:
                f.writelines(html)

            # Copy the css and statics files next to the report
            for f in glob(stylesheets):
                copy(f, self.prefix)
            copytree(os.path.dirname(statics), os.path.join(self.prefix,'statics'), dirs_exist_ok=True)

        else:

            # Open data stream with context manager
            with zipfile.ZipFile(f"{self.prefix}.zip", mode='w') as zip:
                zip.writestr(output, html)
                for stylesheet in glob(stylesheets):
                    zip.write(stylesheet, arcname=os.path.relpath(path=stylesheet, start=ressources))
                for static in glob(statics):
                    zip.write(static, arcname=os.path.relpath(path=static, start=ressources))