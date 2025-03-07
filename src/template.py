from jinja2 import Environment, FileSystemLoader


def render():
    
    env = Environment(loader=FileSystemLoader('src/templates'))

    template = env.get_template("file.html")

    html = template.render()

    with open("report.html",'w') as f:
        f.writelines(html)