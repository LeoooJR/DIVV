from cli import Program
from sys import exit, stderr


def main():

    try:

        Program().launch()

    except Exception as e:

        print(f"Error: {e}", file=stderr)
        exit(1)


if __name__ == "__main__":
    
    main()
