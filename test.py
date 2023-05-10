from sympy import ntt, intt

a = [100, 200, 200, 0, 100, 200, 200, 0]
b = [100, 200, 300, 400]

a_ntt = ntt(a, 3329)
a_intt = intt(a_ntt, 3329)
print(a_ntt)
print(a_intt)
print(ntt(b, 3329))
