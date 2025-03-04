from cli import Programm
from sys import exit


def main():

    prog = Programm()
    try:
        prog.launch()
    except Exception as e:
        print(e)
        exit(1)


if __name__ == "__main__":
    main()
