from cli import EntryPoint
from console import stderr_console
from rich.panel import Panel
from rich import box
from sys import exit


def main():

    try:
        
        EntryPoint().launch()
    # If a SystemExit error occurs, print the error in the standard error stream
    except SystemExit as e:
        # Print raised error in standard error stream
        stderr_console.print(Panel.fit(str(e), box=box.ROUNDED, title="Execution error", subtitle="System exit as 1", highlight=True), style="error")
        # Exit with code 1 by Unix convention
        exit(1)

# Execute only if the script is run directly and not imported
if __name__ == "__main__":
    
    main()
