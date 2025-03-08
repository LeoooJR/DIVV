from jinja2 import Environment, FileSystemLoader

class Report:

    def __init__(self, vcfs: list[str], df):

        self.vcfs = vcfs
        
        self.data = df

    def __str__(self):
        pass

    def create(self):
        
        env = Environment(loader=FileSystemLoader('src/templates'))

        template = env.get_template("file.html")

        html = template.render(vcfs=self.vcfs, df=self.data)

        with open("report.html",'w') as f:
            f.writelines(html)