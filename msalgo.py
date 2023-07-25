import math
from scipy.special import comb
from lfsr import LFSR

#计算1的位数
def bit_count(val):
    cnt = 0
    while val:
        cnt += 1
        val = val & (val - 1)
    return cnt

#汉明距离变换
hamming = {15 : [[]] * 16, 16 : [[]] * 17}
for i in range(1, 1 << 16):
    cnt = bit_count(i)
    hamming[16][cnt].append(i)
    if i < (1 << 15):
        hamming[15][cnt].append(i)


class MS:
    def __init__(self, mask, n, z, p):
        self.mask = mask    #抽头掩码
        self.n = n          #级数
        self.z = list(z)    #密钥流序列
        self.p = p          #相关概率
        self.len_z = len(z) #密钥流长度
        self.t = 0          #抽头数
        for i in range(n):
            if (1 << i) & self.mask:
                self.t += 1
        self.m = round(self.M())        #平均每一位方程数
        self.s_init = self.S(self.t)    #每一位正确的后验概率s

    @staticmethod
    def bit_stream_to_int(a):
        return int(''.join(map(str, a)), 2)

    def M(self):
        return math.log2(self.len_z / (2 * self.n)) * (self.t + 1)

    def S(self, t):
        if t == 1:
            return self.p
        return self.p * self.S(t - 1) + (1 - self.p) * (1 - self.S(t - 1))
    
    #生成校验等式
    def get_eq(self):
        tap = [self.n]
        for i in range(self.n):
            if (self.mask >> i) & 1:
                tap.append(self.n - i - 1)
        tap.reverse()
        eqs = [tap]
        while True:
            if (tap[-1] << 1) >= self.len_z:
                break
            tmp = tap.copy()
            for i, val in enumerate(tmp):
                tmp[i] = val << 1
            eqs.append(tmp)
            tap = tmp
        return eqs

    #算法A得到loc位正确的概率
    def calc_eq(self, eqs, loc):
        shift_eqs = []
        for eq in eqs:
            for pos in eq:
                offset = loc - pos
                if eq[0] + offset < 0 or eq[-1] + offset >= self.len_z:
                    continue
                shift_eqs.append([i + offset for i in eq])
        m = len(shift_eqs)
        if m == 0:
            return 0, 0, 0
        h = 0
        for eq in shift_eqs:
            # print(eq)
            xor_sum = 0
            for i in eq:
                xor_sum ^= self.z[i]
            if xor_sum == 0:
                h += 1
        p1 = comb(m, h) * pow(self.s_init, h) * pow(1 - self.s_init, m - h)
        p0 = comb(m, h) * pow(self.s_init, m - h) * pow(1 - self.s_init, h)
        return m, h, p1 / (p1 + p0)

    #生成矩阵
    def gen_linear_eq(self):
        length = max(self.len_z, self.n)
        tap = []
        for i in range(self.n):
            if (self.mask >> i) & 1:
                tap.append(i + 1)
        eqs = []
        for i in range(self.n):
            eqs.append(1 << i)
        for i in range(self.n, length):
            res = 0
            for j in tap:
                res ^= eqs[i - j]
            eqs.append(res)
        return eqs

    #根据生成矩阵解方程
    @staticmethod
    def solve(assume, n):
        eq_len = len(assume)
        mat = []
        for i in range(eq_len):
            mat.append([0] * n)
        b = [0] * eq_len
        for i in range(eq_len):
            b[i] = assume[i][1]
            for j in range(n):
                mat[i][j] = (assume[i][0] >> j) & 1
        for i in range(n):
            tmp = -1
            for j in range(i, eq_len):
                if mat[j][i]:
                    tmp = j
                    break
            if tmp == -1:
                return []
            mat[tmp], mat[i] = mat[i], mat[tmp]
            b[tmp], b[i] = b[i], b[tmp]
            for j in range(eq_len):
                if not mat[j][i] or i == j:
                    continue
                b[j] ^= b[i]
                for k in range(i, n):
                    mat[j][k] ^= mat[i][k]
        if not any(mat[n - 1]):
            return []
        # print(b[:n])
        return b[:n]

    #计算初始状态并验证
    def get_init_stat(self, locs, linear_eq):
        assume = [(linear_eq[x[0]], x[1]) for x in locs]
        b = []
        idx = self.n
        # print("----- try solve equations -----")
        while not b:
            b = MS.solve(assume[:idx], self.n)
            idx += 1
        # print("----- solve success -----")
        stat = MS.bit_stream_to_int(b)

        # print("----- genrate original LFSR -----")
        l = LFSR(stat, self.mask, self.n)
        for i in range(self.n):
            l.step_back()
        init_stat = l.init
        # print("init:", init_stat)
        # print("----- genrate original LFSR finished -----")

        same_cnt = 0
        for i in range(self.len_z):
            same_cnt += int(self.z[i] == l.next())
        rate = same_cnt / self.len_z
        if abs(rate - self.p) < 0.05:
            return init_stat
        else:
            #hamming变换
            for i in range(self.n + 1):
                for filp in hamming[self.n][i]:
                    change = init_stat ^ filp
                    l.init = change
                    cnt = 0
                    for i in range(self.len_z):
                        cnt += int(self.z[i] == l.next())
                    rate = cnt / self.len_z
                    if abs(rate - self.p) < 0.05:
                        return change
        return init_stat

    #算法A攻击函数
    def crackA(self):
        eqs = self.get_eq()
        # print("----- select candidates -----")
        candidates = []
        for i in range(self.len_z):
            m, h, p_star = self.calc_eq(eqs, i)
            if p_star > 0.5:
                candidates.append((p_star, i, m, h))
        candidates.sort(reverse=True)
        # candidates = candidates[:2*self.n]
        # print(candidates[:5])
        # print("----- select candidates finished -----")
        linear_eq = self.gen_linear_eq()
        locs = [(cand[1], self.z[cand[1]]) for cand in candidates]
        return self.get_init_stat(locs, linear_eq)

    #计算s_init
    @staticmethod
    def var_S(var_p, t):
        assert(t == len(var_p))
        if t == 1:
            return var_p[0]
        s = MS.var_S(var_p[:-1], t - 1)
        return var_p[-1] * s + (1 - var_p[-1]) * (1 - s)

    def Q(self, h):
        res = 0
        for i in range(h + 1):
            res += comb(self.m, i) * (self.p_false(i) + self.p_true(i))
        return res

    def I(self, h):
        res = 0
        for i in range(h + 1):
            res += comb(self.m, i) * (self.p_false(i) - self.p_true(i))
        return res

    def p_true(self, h):
        return self.p * pow(self.s_init, h) * pow(1 - self.s_init, self.m - h)

    def p_false(self, h):
        return (1 - self.p) * pow(1 - self.s_init, h) * pow(self.s_init, self.m - h)

    def p_update(self, h):
        t = self.p_true(h)
        f = self.p_false(h)
        return t / (t + f)

    def calc_h_max(self):
        h_max = 0
        I_max = 0
        for i in range(self.m + 1):
            newI = self.I(i)
            if newI > I_max:
                I_max = newI
                h_max = i
        return h_max

    #计算迭代后验概率的4个部分
    def parcheck(self, eqs, arr_p):
        poly = self.t + 1
        s = [[1, 1, 1, 1] for i in range(self.len_z)]
        for eq in eqs:
            curtaps = eq.copy()
            offset = self.len_z - eq[-1]
            # print('\n' , offset , '\n')
            for i in range(offset):
                # print(i , ' ')
                xor_sum = 0
                for tap in curtaps:
                    xor_sum ^= self.z[tap]
                for j in range(poly):
                    var_p = [arr_p[curtaps[k]] for k in range(poly) if k != j]
                    cur_s = MS.var_S(var_p, self.t)
                    cur_bit = curtaps[j]
                    if xor_sum == 0:
                        s[cur_bit][0] *= cur_s
                        s[cur_bit][1] *= 1 - cur_s
                    else:
                        s[cur_bit][2] *= cur_s
                        s[cur_bit][3] *= 1 - cur_s
                for j in range(poly):
                    curtaps[j] += 1
        return s

    #计算迭代后验概率
    @staticmethod
    def var_p_update(p, arr_s):
        t = p * arr_s[0] * arr_s[3]
        div = t + (1 - p) * arr_s[1] * arr_s[2]
        #问题出在这里
        if div == 0:
            return 1
        return t / (t + (1 - p) * arr_s[1] * arr_s[2])

    #检验是否满足校验方程
    @staticmethod
    def check(eq, z):
        offset = len(z) - eq[-1]
        for i in range((offset)):
            xor_sum = 0
            for tap in eq:
                xor_sum ^= z[tap]
            if xor_sum:
                return False
            for j, val in enumerate(eq):
                eq[j] = val + 1
        return True

    #算法B攻击函数
    def crackB(self):
        prim_z = self.z.copy()
        eqs = self.get_eq()
        h_max = self.calc_h_max()
        p_thr = (self.p_update(h_max) + self.p_update(h_max + 1)) / 2
        N_thr = self.Q(h_max) * self.len_z
        # print("------------------------------------------")
        # print("Round\tIteration\t  N_w")
        # print("------------------------------------------")
        for r in range(1, 1000):
            arr_p = [self.p for i in range(self.len_z)]
            for iter in range(1, 6):
                arr_ss = self.parcheck(eqs, arr_p)
                for i in range(self.len_z):
                    arr_p[i] = MS.var_p_update(arr_p[i], arr_ss[i])
                N_w = 0
                for i in range(self.len_z):
                    if arr_p[i] < p_thr:
                        N_w += 1
                # print(r, '\t', iter, '\t\t', N_w)
                if N_w >= N_thr:
                    break
            cnt_filp_bit = 0
            for pos in range(self.len_z):
                if arr_p[pos] < p_thr:
                    self.z[pos] ^= 1
                    cnt_filp_bit += 1
            # print("------------------------------------------")
            if MS.check(eqs[0][:], self.z) or cnt_filp_bit == 0:
                stat = MS.bit_stream_to_int(self.z[:self.n])
                l = LFSR(stat, self.mask, self.n)
                for i in range(self.n):
                    l.step_back()
                self.z = prim_z
                return l.init
        return list() 