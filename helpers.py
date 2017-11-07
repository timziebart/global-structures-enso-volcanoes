
import traceback
import sys


def printException(text = "", file = sys.stderr):

    print(file = file)
    print("#"*40, file = file)
    exc_type, exc_value, exc_traceback = sys.exc_info()
    lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
    print(''.join('!! ' + line for line in lines) + "\n", file = sys.stderr)

    if text:
        print(text + "\n", file = sys.stderr)
    print("#"*40, file = file)
    print(file = file)