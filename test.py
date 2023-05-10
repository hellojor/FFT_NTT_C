import numpy as np

a = np.array([1, 2, 3, 4])
a *= 2
print(np.fft.ifft(a))
