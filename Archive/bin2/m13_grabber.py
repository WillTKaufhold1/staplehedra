def m13seq():
    from functools import reduce
    with open('../m13mp18.dat','r') as f:
        data = f.readlines()
    filtered = list(filter(lambda x : x[0] != ';',data))
    stripped = list(map(lambda x : x.replace(' ','').replace('\n',''), filtered))
    accumulalted = str(reduce( lambda x,y : x + y, stripped, ''))

    return accumulalted