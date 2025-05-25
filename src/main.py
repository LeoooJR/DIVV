from cli import Program
from console import stderr_console
from rich.panel import Panel
from rich import box
from sys import exit


def main():

    try:

        Program().launch()

    except SystemExit as e:

        stderr_console.print(Panel.fit(e, box=box.ROUNDED, title="Execution error", subtitle="System exit as 1"), style="error")
        exit(1)

if __name__ == "__main__":
    
    main()
