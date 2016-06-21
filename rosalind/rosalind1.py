

if __name__ == "__main__":
    print 'Holla World'

    a = 851
    b = 894

    output = a**2 + b**2

    print output

    fn = 'camelot.txt'
    with open(fn, 'r') as f:
        for anything in f:
            print anything
