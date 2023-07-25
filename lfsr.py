import random

class LFSR:
    def __init__(self, init, mask, length):
        self.init = init
        self.length = length
        self.lengthmask = (1 << length) - 1
        self.mask = mask & self.lengthmask

    def next(self):
        nextdata = (self.init << 1) & self.lengthmask
        output = LFSR.parity(self.init & self.mask)
        nextdata ^= output
        self.init = nextdata
        return output

    def step_back(self):
        output = self.init & 1
        predata = self.init >> 1
        high_bit = LFSR.parity(predata & self.mask) ^ output
        self.init = (high_bit << (self.length - 1)) | predata

    @staticmethod
    def parity(x):
        res = 0
        while x:
            x -= x & (-x)
            res ^= 1
        return res

def combine(x1,x2,x3):
    return (x1*x2)^(x2*x3)^(x1*x3)

if __name__=="__main__":
    # length = 48
    # mask = (1 << length)
    # seed = [random.randint(0, (1 << length) - 1) for i in range(3)]
    # print(seed)
    # l1, l2, l3 = [lfsr(seed[i],0b100000000000000000000000010000000000000000000000,48) for i in range(3)]
    l1 = lfsr(173975563016619,0b100000000000000000000000010000000000000000000000,48)
    l2 = lfsr(170456896330817,0b100000000000000000000000000000000010000000000000,48)
    l3 = lfsr(254656305199092,0b100000100000000000000000000000000000000000000000,48)
    txt = []
    with open("keystream","wb") as f:
        for i in range(8192):
            b = 0
            tmp = 0
            for j in range(8):
                tmp = combine(l1.next(),l2.next(),l3.next())
                b = (b<<1)+tmp
                txt.append(str(tmp))
            f.write(chr(b).encode())
    with open("keystream.txt","w") as f:
        f.write(''.join(txt))