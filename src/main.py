from cli import Program
from sys import exit, stderr


def main():

    prog = Program()
    try:
        prog.launch()
    except Exception as e:
        print(f"Error: {e.args[0]}", file=stderr)
        exit(1)


if __name__ == "__main__":
    main()
