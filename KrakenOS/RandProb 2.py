
import random

def prov(pro):
    """prov.

    Parameters
    ----------
    pro :
        pro
    """
    a_list = [0, 1]
    prob = pro
    distribution = [prob, (1.0 - prob)]
    random_number = random.choices(a_list, distribution)
    return random_number
for i in range(0, 100):
    rnum = prov(0.05)
    print(rnum)

