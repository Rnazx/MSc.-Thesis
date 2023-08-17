import os
print(os.getcwd())
def f(x):
    return x*x

def fn():    
    from multiprocessing import Pool
    with Pool(5) as p:
        return p.map(f, [1, 2, 3])
if __name__ == "__main__":
    print(fn())